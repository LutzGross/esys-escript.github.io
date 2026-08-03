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
Checks that toFinley() carries the tag NAMES and the Dirac points across.

Neither travels with the mesh arrays: tag names are domain state, and a Dirac
point is resolved by oxley's addPoints() to a node rather than being a mesh
entity. Both used to be dropped silently - the export simply had four
hardcoded boundary names and an empty Point table.

The Dirac check is the interesting one. oxley snaps each point to the nearest
OWNED lnodes node, so the exported mesh must place it at THAT node, not at the
coordinates the user asked for and not at whatever node of the export happens
to be nearest - the export has more nodes than oxley (it materialises the
hanging positions), so a fresh search could pick one oxley would never choose.

Positions are compared by moments summed over ranks, because which rank owns a
point may legitimately differ between the two domains: oxley assigns it to the
owner of the nearest node, finley to the owner of the DOF after prepare().

Run with: run-escript -n<P> -t1 validate_export_tags_dirac.py
"""

import sys

from esys.escript import *
from esys.oxley import Rectangle
import esys.finley

RANK = getMPIRankWorld()

# asked-for points, deliberately NOT on nodes, so that snapping is visible
POINTS = [(0.3, 0.3), (0.7, 0.7), (0.1, 0.9)]
TAGS = ["A", "B", "C"]

CASES = [
    ("uniform level 2", 2),
    ("one seam",        [[2], [3]]),
    ("3x3 mixed",       [[3, 1, 2], [1, 2, 1], [2, 1, 3]]),
]

failures = []


def check(ok, what):
    if not ok:
        failures.append(what)
    if RANK == 0:
        print("    %-58s %s" % (what, "OK" if ok else "*** FAILED ***"))
        sys.stdout.flush()


def moments(fs):
    """count and coordinate moments of a Dirac function space, summed globally.

    A fingerprint of the point SET that does not depend on which rank holds
    which point, so the two domains can be compared without gathering.

    Fixed point, because getMPIWorldSum only reduces ints. The two domains hold
    the same doubles - finley copies the coordinates oxley resolved - so a 1e-6
    grid is a strict comparison here, not a tolerance.
    """
    x = fs.getX()
    n = x.getNumberOfDataPoints()
    s = [0] * 4
    for i in range(n):
        p = x.getTupleForDataPoint(i)
        s[0] += int(round(p[0] * 1e6))
        s[1] += int(round(p[1] * 1e6))
        s[2] += int(round((p[0] * p[0] + p[1] * p[1]) * 1e6))
        s[3] += int(round(p[0] * p[1] * 1e6))
    # every one of these is collective, so every rank must reach all of them
    return [getMPIWorldSum(n)] + [getMPIWorldSum(v) for v in s]


def hanging_position(dom):
    """
    One hanging position of the forest, agreed by every rank, or None.

    This is the adversarial place to ask for a Dirac point. It is NOT a node of
    the oxley domain, so oxley snaps to an lnodes node half an element away -
    but it IS a node of the export, at distance zero. A converter that located
    points on the exported mesh instead of reusing oxley's answer would put the
    point here and disagree with oxley, and nothing else in this script would
    notice.

    Picked as the lexicographic maximum over ranks so that all agree without a
    broadcast: getMPIWorldMax only reduces ints, hence fixed point again, and
    the y is taken only from ranks that achieved the winning x so the two
    coordinates cannot come from different points.
    """
    info = dom.getMeshInfo(True)
    xy = info["nodeCoords"].reshape(-1, 2)
    best = (-1, -1)
    for i in info["constrainedNodes"]:
        p = xy[int(i)]
        best = max(best, (int(round(p[0] * 1e6)), int(round(p[1] * 1e6))))
    gx = getMPIWorldMax(best[0])
    gy = getMPIWorldMax(best[1] if best[0] == gx else -1)
    if gx < 0:
        return None
    return (gx * 1e-6, gy * 1e-6)


def report(name, levels):
    if RANK == 0:
        print("\n%s" % name)
    if isinstance(levels, int):
        n0 = n1 = 2
    else:
        n0, n1 = len(levels), len(levels[0])

    points, tags = list(POINTS), list(TAGS)
    probe = Rectangle(n0=n0, n1=n1, l0=1., l1=1., refine_level=levels)
    hp = hanging_position(probe)            # collective
    if hp is not None:
        points.append(hp)
        tags.append("H")
        if RANK == 0:
            print("    (asking for a point ON the hanging position %s)" % (hp,))

    dom = Rectangle(n0=n0, n1=n1, l0=1., l1=1., refine_level=levels,
                    diracPoints=points, diracTags=tags)
    dom.setTagMap("myregion", 7)
    fin = dom.toFinley()

    # --- tag names -------------------------------------------------------
    names = [s.strip() for s in dom.showTagNames().split(",") if s.strip()]
    missing = [t for t in names if not fin.isValidTagName(t)]
    check(not missing, "every oxley tag name survives (missing: %s)" % missing)
    wrong = [t for t in names
             if fin.isValidTagName(t) and fin.getTag(t) != dom.getTag(t)]
    check(not wrong, "and keeps its value (differing: %s)" % wrong)
    check(fin.isValidTagName("myregion") and fin.getTag("myregion") == 7,
          "a user-defined name is carried, not just the boundary ones")

    # --- Dirac points ----------------------------------------------------
    mo = moments(DiracDeltaFunctions(dom))
    mf = moments(DiracDeltaFunctions(fin))
    check(mf[0] == len(points),
          "export has all %d Dirac points (has %d)" % (len(points), mf[0]))
    check(mo[0] == mf[0],
          "point count matches oxley (%d vs %d)" % (mo[0], mf[0]))
    same = (mo[1:] == mf[1:])
    check(same, "points sit exactly where oxley put them")

    # a Dirac point must never land on a materialised hanging position: those
    # are not nodes of the oxley domain, so oxley could not have chosen one
    for tag in tags:
        d = Data(0., DiracDeltaFunctions(fin))
        d.setTaggedValue(tag, 5.)
        peak = sup(d)                      # collective
        check(peak == 5., "tag '%s' selects a point on the export" % tag)


if RANK == 0:
    print("=== tag names and Dirac points through toFinley(), %d rank(s) ==="
          % getMPISizeWorld())

for name, levels in CASES:
    report(name, levels)

if RANK == 0:
    print("\nRESULT: %s" % ("ALL CHECKS PASSED" if not failures
                            else "FAILED (%s)" % ", ".join(failures)))
