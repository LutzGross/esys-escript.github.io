##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################
"""
Validation of the oxley VTK/Silo output (design milestone A3).

Writes a uniform, conforming Block domain via the escript-standard weipa path
(saveVTK / saveSilo) and checks the .vtu structurally against the analytic mesh:
point count, cell count, cell type (VTK_QUAD / VTK_HEXAHEDRON) and that the
written point coordinates match the regular grid.
"""
import glob
import os
import sys
import numpy as np
import xml.etree.ElementTree as ET

from esys.oxley import Block
from esys.escript import ContinuousFunction
import esys.escript as esc
from esys.weipa import saveVTK, saveSilo

WORK = os.environ.get("OXLEY_WORKDIR", ".")
VTK_QUAD, VTK_HEX = 9, 12


def _da(piece, name):
    for da in piece.iter("DataArray"):
        if da.get("Name") == name:
            return np.fromstring(da.text, sep=" ")
    return None


def parse_vtu(path):
    piece = ET.parse(path).getroot().find(".//Piece")
    npts = int(piece.get("NumberOfPoints"))
    ncells = int(piece.get("NumberOfCells"))
    types = _da(piece, "types")
    pts = None
    for pts_el in piece.iter("Points"):
        pts = _da(pts_el, None) if _da(pts_el, None) is not None else None
        for da in pts_el.iter("DataArray"):
            pts = np.fromstring(da.text, sep=" ").reshape(-1, 3)
    return npts, ncells, types, pts


def _check(ok, name, cond):
    print(("  PASS  " if cond else "  FAIL  ") + name)
    ok.append(bool(cond))


def validate(numBlocks, length, origin, refine_level):
    dim = len(numBlocks)
    N = [numBlocks[i] * (2 ** refine_level) for i in range(dim)]
    exp_nnode = int(np.prod([n + 1 for n in N]))
    exp_nelem = int(np.prod(N))
    exp_type = VTK_QUAD if dim == 2 else VTK_HEX
    print("Block numBlocks=%s refine_level=%d -> %d nodes, %d cells (%s)" %
          (numBlocks, refine_level, exp_nnode, exp_nelem,
           "quad" if dim == 2 else "hex"))

    dom = Block(numBlocks=numBlocks, length=length, origin=origin,
                refine_level=refine_level)
    x = ContinuousFunction(dom).getX()

    base = os.path.join(WORK, "oxley_a3_%dd_L%d" % (dim, refine_level))
    # One rank clears the old output, then everyone waits. Removing on every
    # rank is a race - the others find the file already gone and the run dies
    # with FileNotFoundError - and writing before the removal has finished would
    # delete the file just written.
    if esc.getMPIRankWorld() == 0:
        for f in glob.glob(base + "*"):
            os.remove(f)
    esc.MPIBarrierWorld()
    saveVTK(base + ".vtu", u=x)
    saveSilo(base + ".silo", u=x)
    esc.MPIBarrierWorld()               # the checks below read what was written

    ok = []
    vtu = glob.glob(base + "*.vtu")
    _check(ok, "saveVTK produced a .vtu file", len(vtu) == 1)
    silo = glob.glob(base + "*.silo")
    _check(ok, "saveSilo produced a non-empty .silo file",
           len(silo) == 1 and os.path.getsize(silo[0]) > 0)
    if not vtu:
        print("  => FAILED\n")
        return False

    npts, ncells, types, pts = parse_vtu(vtu[0])
    _check(ok, "NumberOfPoints == %d" % exp_nnode, npts == exp_nnode)
    _check(ok, "NumberOfCells == %d" % exp_nelem, ncells == exp_nelem)
    _check(ok, "all cell types == %s" % ("VTK_QUAD" if dim == 2 else "VTK_HEX"),
           types is not None and len(types) == exp_nelem and np.all(types == exp_type))
    if pts is not None:
        grid_lo = np.array(origin)
        grid_hi = np.array(origin) + np.array(length)
        got_lo = pts[:, :dim].min(axis=0)
        got_hi = pts[:, :dim].max(axis=0)
        _check(ok, "point coordinate span matches the box",
               np.allclose(got_lo, grid_lo) and np.allclose(got_hi, grid_hi))
        h = [length[i] / N[i] for i in range(dim)]
        ongrid = all(np.allclose(np.round((pts[:, i] - origin[i]) / h[i]),
                                 (pts[:, i] - origin[i]) / h[i]) for i in range(dim))
        _check(ok, "all points lie on the regular grid", ongrid)

    # data correctness: the written field u (== getX) must equal the point
    # coordinates, i.e. data samples map to the correct nodes.
    u = _da(ET.parse(vtu[0]).getroot().find(".//Piece"), "u")
    if u is not None and pts is not None:
        u = u.reshape(-1, 3)
        _check(ok, "written field u maps to the correct nodes (u == coords)",
               np.allclose(u[:, :dim], pts[:, :dim]))

    print("  => %s\n" % ("OK" if all(ok) else "FAILED"))
    return all(ok)


CASES = [
    ((2, 2),    (1., 1.),      (0., 0.),      0),
    ((2, 2),    (1., 1.),      (0., 0.),      2),
    ((3, 4),    (6., 8.),      (-1., -2.),    1),
    ((1, 1),    (2., 3.),      (5., 5.),      3),
    ((2, 2, 2), (1., 1., 1.),  (0., 0., 0.),  1),
    ((3, 2, 2), (3., 2., 2.),  (-1., 0., -5.), 1),
    ((1, 1, 1), (1., 1., 1.),  (0., 0., 0.),  2),
]

if __name__ == "__main__":
    allok = all(validate(*c) for c in CASES)
    print("ALL PASSED" if allok else "SOME CASES FAILED")
    sys.exit(0 if allok else 1)
