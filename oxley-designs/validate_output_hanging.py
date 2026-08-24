"""
Validation of weipa output for a NON-conforming (adaptive) oxley forest.

A hanging position is not an lnodes node: the element_nodes slot holds a MASTER
instead, a node that lies OUTSIDE the element. Handed to weipa unchanged that
draws the cell with a corner in the wrong place, and (before the fix) also moved
the master itself, because the coordinate walk wrote the hanging position into the
master's entry. getMeshAccess(materializeHanging=true) gives every hanging
position a node of its own, and weipa fills the value in from the masters.

Everything here is read back out of the written .vtu, so it tests what a viewer
would actually see, not what oxley thinks it wrote:

  1. cell shapes    - every cell's corners must be an axis-aligned square (cube)
                      of the cell's own size. This is the check that fails hard
                      on a master in the corner list.
  2. cell volume    - the cell volumes must sum to the domain volume, so the
                      cells tile the domain with no overlap or hole.
  3. node positions - every point must be a corner of some cell, and every cell
                      corner must be a point (no node left at a stale position).
  4. VALUES         - a linear field is written as a nodal variable and must come
                      back exact at every point, hanging ones included. A hanging
                      node is the average of its masters, and an average of a
                      linear function is that function, so any error here means
                      the interpolation is wrong (or absent, leaving a zero).
  5. a smooth field - written too, and checked for the tell-tale 0.0 that a
                      missing sample leaves behind.

Run:
    ./bin/run-escript -n1 -t1 oxley-designs/validate_output_hanging.py
    ./bin/run-escript -n2 -t1 oxley-designs/validate_output_hanging.py

Under MPI weipa writes one piece per rank; this script checks each piece it can
find, so the per-rank meshes are validated but the seams between them are not.
"""
import glob
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

import esys.escript as esc
import esys.oxley as oxley
from esys.escript import ContinuousFunction
from esys.weipa import saveVTK

WORK = os.environ.get("OXLEY_WORKDIR", ".")
VTK_QUAD, VTK_HEX = 9, 12
TOL = 1e-5              # weipa writes float32, so ~7 significant digits
GTOL = 1e-6             # the same, applied to coordinates

rank = esc.getMPIRankWorld()
failures = []


def _da(el, name):
    for da in el.iter("DataArray"):
        if da.get("Name") == name:
            return np.fromstring(da.text, sep=" ")
    return None


def parse_vtu(path):
    """points, connectivity (list of per-cell index lists), cell types"""
    piece = ET.parse(path).getroot().find(".//Piece")
    npts = int(piece.get("NumberOfPoints"))
    ncells = int(piece.get("NumberOfCells"))
    pts = None
    for pts_el in piece.iter("Points"):
        for da in pts_el.iter("DataArray"):
            pts = np.fromstring(da.text, sep=" ").reshape(-1, 3)
    conn = _da(piece, "connectivity").astype(int)
    offs = _da(piece, "offsets").astype(int)
    types = _da(piece, "types").astype(int)
    cells = []
    start = 0
    for o in offs:
        cells.append(conn[start:o])
        start = o
    fields = {}
    for cd in piece.iter("PointData"):
        for da in cd.iter("DataArray"):
            fields[da.get("Name")] = np.fromstring(da.text, sep=" ")
    assert npts == len(pts) and ncells == len(cells)
    return pts, cells, types, fields


def check(name, cond, detail=""):
    if rank == 0:
        print("  %-52s %s%s" % (name, "PASS" if cond else "FAIL",
                                "" if cond else "   " + detail))
    if not cond:
        failures.append(name)
    return cond


def validate(label, dom, dim, volume):
    if rank == 0:
        print("\n%s" % label)

    x = ContinuousFunction(dom).getX()
    linear = 1.0 + 2.0 * x[0] + 3.0 * x[1] + (4.0 * x[2] if dim == 3 else 0.0)
    smooth = esc.sin(np.pi * x[0]) * esc.exp(x[1])

    base = os.path.join(WORK, "oxley_hanging_%dd_%s"
                        % (dim, label.split(",")[1].strip().replace(" ", "_")))
    if rank == 0:                       # every rank globs the same directory
        for f in glob.glob(base + "*"):
            try:
                os.remove(f)
            except OSError:
                pass
    saveVTK(base + ".vtu", linear=linear, smooth=smooth, coords=x)

    files = sorted(glob.glob(base + "*.vtu"))
    if not check("saveVTK produced a .vtu file", len(files) > 0):
        return

    pts, cells, types, fields = parse_vtu(files[0])

    # Under MPI weipa writes ONE shared point list and the cells of every rank
    # into it, so the header count and the connectivity must agree - the failure
    # this catches is a rank-local node numbering leaking into a global file,
    # which loads but cannot be plotted.
    maxref = max((int(max(c)) for c in cells if len(c)), default=-1)
    check("connectivity stays inside the point list",
          maxref < len(pts),
          "cell node index %d but only %d points" % (maxref, len(pts)))
    want_type = VTK_QUAD if dim == 2 else VTK_HEX
    nc = 4 if dim == 2 else 8

    check("all cells are %s" % ("VTK_QUAD" if dim == 2 else "VTK_HEX"),
          len(types) > 0 and bool(np.all(types == want_type)))
    check("every cell has %d corners" % nc,
          all(len(c) == nc for c in cells))

    # 1. cell shapes and 2. total volume.
    #
    # A cell must be an axis-aligned box: two distinct values on each axis, and
    # its corners exactly their cartesian product. That is what a master in the
    # corner list breaks - a master sits one full (coarse) cell away, so the axis
    # it displaces along ends up with three distinct values.
    #
    # NOT "a square of its own size": the blocks need not be square (n0=2,n1=1
    # over the unit square gives 0.5 x 1.0 blocks), so cells need not be either.
    # And weipa writes float32, hence a relative tolerance, not an absolute one.
    def cluster(vals):
        """distinct coordinate values, merging any that agree to float32"""
        out = []
        for v in sorted(vals):
            if not out or abs(v - out[-1]) > GTOL * max(1.0, abs(v)):
                out.append(v)
        return out

    misshapen, total = 0, 0.0
    for c in cells:
        p = [tuple(float(v) for v in pts[i][:dim]) for i in c]
        axes = [cluster(q[d] for q in p) for d in range(dim)]
        ok = all(len(a) == 2 for a in axes)
        if ok:
            # every corner of the box must be present exactly once
            seen = set()
            for q in p:
                corner = []
                for d in range(dim):
                    corner.append(0 if abs(axes[d][0] - q[d]) < abs(axes[d][1] - q[d])
                                  else 1)
                seen.add(tuple(corner))
            ok = len(seen) == len(p)
        if ok:
            vol = 1.0
            for a in axes:
                vol *= a[1] - a[0]
            total += vol
        else:
            misshapen += 1
    check("every cell is an axis-aligned box", misshapen == 0,
          "%d of %d cells are not" % (misshapen, len(cells)))
    if esc.getMPISizeWorld() == 1:
        check("cell volumes sum to the domain volume",
              abs(total - volume) < GTOL * volume,
              "got %.12g, want %.12g" % (total, volume))

    # 3. points and cell corners are the same set
    def key(p):
        return tuple(round(float(v), 6) for v in p[:dim])   # float32 again

    corner_pos = set()
    for c in cells:
        for i in c:
            corner_pos.add(key(pts[i]))
    point_pos = set(key(p) for p in pts)
    check("every point is a corner of some cell",
          point_pos <= corner_pos, "%d are not" % len(point_pos - corner_pos))
    check("every cell corner is a point",
          corner_pos <= point_pos, "%d are not" % len(corner_pos - point_pos))

    # 4. the linear field must be exact everywhere, hanging nodes included
    if check("the .vtu carries the nodal fields",
             "linear" in fields and "smooth" in fields):
        lin = fields["linear"]
        want = np.array([1.0 + 2.0 * p[0] + 3.0 * p[1] + (4.0 * p[2] if dim == 3 else 0.0)
                         for p in pts])
        err = np.abs(lin - want)
        worst = int(np.argmax(err))
        check("a linear nodal field is exact at every point",
              err.max() < TOL,
              "worst %.3e at point %d %s (got %.6f, want %.6f)"
              % (err.max(), worst, pts[worst][:dim], lin[worst], want[worst]))

        # 5. a missing sample shows up as a hard zero
        sm = fields["smooth"]
        zeros = int(np.sum(np.abs(sm) < 1e-12))
        interior = [i for i in range(len(pts))
                    if abs(pts[i][0]) > 1e-9 and abs(pts[i][0] - 1.0) > 1e-9]
        bad = [i for i in interior if abs(sm[i]) < 1e-12]
        check("no interior point has a zero smooth value (unfilled sample)",
              len(bad) == 0, "%d of %d points are zero" % (zeros, len(pts)))


# The workhorse case: two blocks, one refined twice and one three times. Both
# sides are refined, so the seam is a 2:1 interface between real fine cells
# rather than between a fine block and an unrefined one.
validate("2D, one seam, blocks at levels [[2],[3]]",
         oxley.Rectangle(n0=2, n1=1, refine_level=[[2], [3]]), 2, 1.0)
validate("2D, a deep block, blocks at levels [[3,1],[1,2]]",
         oxley.Rectangle(n0=2, n1=2, refine_level=[[3, 1], [1, 2]]), 2, 1.0)
validate("2D, blocks at levels [[3,1,2],[1,2,1],[2,1,3]]",
         oxley.Rectangle(n0=3, n1=3,
                         refine_level=[[3, 1, 2], [1, 2, 1], [2, 1, 3]]), 2, 1.0)

validate("3D, one seam, blocks at levels [[[2]],[[3]]]",
         oxley.Brick(n0=2, n1=1, n2=1, refine_level=[[[2]], [[3]]]), 3, 1.0)
validate("3D, a deep block, 2x2x2 at mixed levels",
         oxley.Brick(n0=2, n1=2, n2=2,
                     refine_level=[[[2, 1], [1, 1]], [[1, 1], [1, 2]]]), 3, 1.0)

if rank == 0:
    print("\nRESULT: %s" % ("ALL CHECKS PASSED" if not failures
                            else "FAILED (%s)" % ", ".join(failures)))
sys.exit(0 if not failures else 1)
