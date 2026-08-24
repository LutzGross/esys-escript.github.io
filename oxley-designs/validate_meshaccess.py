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
Validation of the oxley lnodes mesh-access interface (design milestone A2).

For a uniform, conforming Block domain the mesh is a regular grid, so the
lnodes-based getMeshInfo() output can be checked analytically:

  N_i          = numBlocks_i * 2**refine_level        (elements per axis)
  numElements  = prod(N_i)
  numNodes     = prod(N_i + 1)                         (shared nodes deduplicated)
  nodeCoords   = the regular grid points
  connectivity = a consistent coordinate<->index bijection

The node-count / dedup checks are what distinguish correct lnodes numbering from
the old floating-point coordinate hashing.
"""
import itertools
import sys
import numpy as np
from esys.oxley import Block


def _check(ok_list, name, cond):
    print(("  PASS  " if cond else "  FAIL  ") + name)
    ok_list.append(bool(cond))


def _coord_set(a):
    return set(tuple(r) for r in np.round(a, 9))


def validate(numBlocks, length, origin, refine_level):
    dim = len(numBlocks)
    N = [numBlocks[i] * (2 ** refine_level) for i in range(dim)]
    h = [length[i] / N[i] for i in range(dim)]
    nc = 4 if dim == 2 else 8
    exp_nelem = int(np.prod(N))
    exp_nnode = int(np.prod([n + 1 for n in N]))
    print("Block numBlocks=%s length=%s origin=%s refine_level=%d -> %s grid, "
          "%d elements, %d nodes" % (numBlocks, length, origin, refine_level,
          "x".join(map(str, N)), exp_nelem, exp_nnode))

    dom = Block(numBlocks=numBlocks, length=length, origin=origin,
                refine_level=refine_level)
    info = dom.getMeshInfo()
    coords = info["nodeCoords"]
    conn = info["elementNodes"]

    ok = []
    _check(ok, "numDim", info["numDim"] == dim)
    _check(ok, "nodesPerElement == %d" % nc, info["nodesPerElement"] == nc)
    _check(ok, "numElements == prod(N) = %d" % exp_nelem,
           info["numElements"] == exp_nelem)
    _check(ok, "numNodes == prod(N+1) = %d  (dedup)" % exp_nnode,
           info["numNodes"] == exp_nnode)
    _check(ok, "nodeCoords shape (%d,%d)" % (exp_nnode, dim),
           coords.shape == (exp_nnode, dim))
    _check(ok, "elementNodes shape (%d,%d)" % (exp_nelem, nc),
           conn.shape == (exp_nelem, nc))
    _check(ok, "serial global ids are identity",
           np.array_equal(info["nodeLnodesId"], np.arange(exp_nnode)))
    _check(ok, "connectivity indices in [0,numNodes)",
           conn.size and conn.min() >= 0 and conn.max() < exp_nnode)
    _check(ok, "every node referenced by an element",
           len(np.unique(conn)) == exp_nnode)

    for d in range(dim):
        lo, hi = float(coords[:, d].min()), float(coords[:, d].max())
        _check(ok, "axis %d span == [%g,%g]" % (d, origin[d], origin[d] + length[d]),
               np.isclose(lo, origin[d]) and np.isclose(hi, origin[d] + length[d]))

    axes = [origin[d] + h[d] * np.arange(N[d] + 1) for d in range(dim)]
    grid = np.array(list(itertools.product(*axes)))
    _check(ok, "node coordinates == regular grid", _coord_set(coords) == _coord_set(grid))
    _check(ok, "no duplicate node coordinates",
           np.unique(np.round(coords, 9), axis=0).shape[0] == exp_nnode)

    cells = coords[conn]                       # (nelem, nc, dim)
    sizes = cells.max(axis=1) - cells.min(axis=1)
    _check(ok, "every element is an axis-aligned cell of size h",
           np.allclose(sizes, h))

    _check(ok, "element tags default to 0",
           np.array_equal(info["elementTags"], np.zeros(exp_nelem, dtype=info["elementTags"].dtype)))

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
