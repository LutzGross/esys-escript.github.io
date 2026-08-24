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
The 3D mesh view of a GRADED forest: does it name every hanging position, and
do the ranks agree on what to call each one?

This is the layer under the tetrahedral split. The split cannot be right if the
view is wrong, and the two failures look nothing alike: a wrong split shows up
as a cracked mesh, a wrong view as a node in the wrong place or two ranks
disagreeing about a node's identity, which under MPI is a crack finley cannot
even see.

What is checked, per forest and per rank:

  BOX          every element's eight node coordinates form the axis-aligned box
               the octant is, corner c at the z-order combination of the two
               extents. A hanging slot that was not redirected still holds a
               MASTER, which lies outside the octant - so this one check catches
               the whole class.
  AVERAGE      a materialised node sits at the average of its masters.
  TABLES       elementFaceHangingNode and elementEdgeHangingNode name nodes at
               the centre of that face and the midpoint of that edge.
  COMPLETE     the converse: wherever a node exists at a face centre or an edge
               midpoint of a local element, the tables name it. This is what
               catches a seam the iterator never reported - the coarse side
               would simply not know the node is there and would leave a hole
               in the export.
  IDS          the export ids are unique on a rank, and every rank that holds a
               given position derives the SAME id for it. The whole numbering
               exists to make that true without an agreement protocol.
"""

import sys
import numpy as np

from esys.escript import getMPISizeWorld, getMPIRankWorld
from esys.oxley import Brick

sys.path.insert(0, "oxley/test/python")
import oxley_meshes

TOL = 1e-9

# p8est's own tables, as the mesh view uses them
EDGE_CORNERS = [(0,1),(2,3),(4,5),(6,7),(0,2),(1,3),(4,6),(5,7),
                (0,4),(1,5),(2,6),(3,7)]
# face f fixes axis f//2 at the low (f even) or high (f odd) extent
FACE_AXIS = [(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)]


def key(p):
    """a position quantised so it can be a dictionary key"""
    return tuple(int(round(c * 1e9)) for c in p)


def check(name, levels):
    rank = getMPIRankWorld()
    dom = Brick(n0=len(levels), n1=len(levels[0]), n2=len(levels[0][0]),
                l0=1., l1=1., l2=1., refine_level=levels)
    info = dom.getMeshInfo(True)

    X = info["nodeCoords"].reshape(-1, 3)
    EN = info["elementNodes"]
    ne = int(info["numElements"])
    nreal = int(info["numRealNodes"])
    fid = info["nodeFinleyId"]
    fh = info["elementFaceHangingNode"].reshape(ne, 6)
    eh = info["elementEdgeHangingNode"].reshape(ne, 12)
    cn = info["constrainedNodes"]
    cm = info["constraintMasters"]
    cw = info["constraintWeights"]

    ok = []

    # ---- BOX ---------------------------------------------------------------
    bad = 0
    lo = np.zeros((ne, 3))
    hi = np.zeros((ne, 3))
    for e in range(ne):
        P = X[EN[e]]
        lo[e] = P.min(axis=0)
        hi[e] = P.max(axis=0)
        for c in range(8):
            want = [lo[e][d] if not ((c >> d) & 1) else hi[e][d]
                    for d in range(3)]
            if np.abs(P[c] - want).max() > TOL:
                bad += 1
    ok.append(("every octant's eight corners form its box", bad == 0, bad))

    # ---- AVERAGE -----------------------------------------------------------
    bad = 0
    mpc = int(info["mastersPerConstrainedNode"])
    for i, node in enumerate(cn):
        w = cw[i]
        mids = cm[i]
        p = np.zeros(3)
        for k in range(mpc):
            if mids[k] >= 0:
                p += w[k] * X[mids[k]]
        if np.abs(X[node] - p).max() > TOL:
            bad += 1
    ok.append(("a materialised node is the average of its masters",
               bad == 0, bad))

    # ---- TABLES ------------------------------------------------------------
    badf = bade = 0
    for e in range(ne):
        centre = 0.5 * (lo[e] + hi[e])
        for f in range(6):
            if fh[e][f] < 0:
                continue
            want = centre.copy()
            axis, side = FACE_AXIS[f]
            want[axis] = hi[e][axis] if side else lo[e][axis]
            if np.abs(X[fh[e][f]] - want).max() > TOL:
                badf += 1
        for ed in range(12):
            if eh[e][ed] < 0:
                continue
            a, b = EDGE_CORNERS[ed]
            want = 0.5 * (X[EN[e][a]] + X[EN[e][b]])
            if np.abs(X[eh[e][ed]] - want).max() > TOL:
                bade += 1
    ok.append(("elementFaceHangingNode sits at the face centre", badf == 0, badf))
    ok.append(("elementEdgeHangingNode sits at the edge midpoint", bade == 0, bade))

    # ---- COMPLETE ----------------------------------------------------------
    # every node this rank holds, by position
    at = {}
    for i in range(len(X)):
        at[key(X[i])] = i
    missf = misse = 0
    for e in range(ne):
        centre = 0.5 * (lo[e] + hi[e])
        for f in range(6):
            want = centre.copy()
            axis, side = FACE_AXIS[f]
            want[axis] = hi[e][axis] if side else lo[e][axis]
            there = at.get(key(want), -1)
            if there >= 0 and fh[e][f] != there:
                missf += 1
        for ed in range(12):
            a, b = EDGE_CORNERS[ed]
            want = 0.5 * (X[EN[e][a]] + X[EN[e][b]])
            there = at.get(key(want), -1)
            if there >= 0 and eh[e][ed] != there:
                misse += 1
    ok.append(("no face centre is left out of the table", missf == 0, missf))
    ok.append(("no edge midpoint is left out of the table", misse == 0, misse))

    # ---- IDS ---------------------------------------------------------------
    ids = [int(v) for v in fid]
    ok.append(("the export ids are unique on this rank",
               len(set(ids)) == len(ids), len(ids) - len(set(ids))))

    disagree = 0
    if getMPISizeWorld() > 1:
        from mpi4py import MPI
        mine = {key(X[i]): int(fid[i]) for i in range(len(X))}
        for other in MPI.COMM_WORLD.allgather(mine):
            for k, v in other.items():
                if k in mine and mine[k] != v:
                    disagree += 1
    ok.append(("every rank derives the same id for a shared position",
               disagree == 0, disagree))

    allok = all(o[1] for o in ok)
    if rank == 0:
        print("%-10s %d elements, %d nodes (%d materialised)"
              % (name, ne, len(X), len(X) - nreal))
        for what, good, count in ok:
            print("  %-4s %-55s %s" % ("PASS" if good else "FAIL", what,
                                       "" if good else "(%d bad)" % count))
    return allok


CASES = [("mixed", oxley_meshes.MIXED_3D),
         ("seam", oxley_meshes.SEAM_3D),
         ("isolated", oxley_meshes.ISOLATED_3D),
         ("peak", oxley_meshes.PEAK_3D),
         ("uniform", [[[1, 1], [1, 1]], [[1, 1], [1, 1]]])]

if __name__ == "__main__":
    good = all([check(n, l) for n, l in CASES])
    if getMPIRankWorld() == 0:
        print("\n" + ("ALL PASSED" if good else "SOME CHECKS FAILED"))
    sys.exit(0 if good else 1)
