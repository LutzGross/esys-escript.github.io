"""
2D oxley -> finley export of a NON-conforming (adaptive) forest.

finley cannot represent a hanging node: one ReferenceElementSet per ElementFile,
and escript's q/r is pointwise Dirichlet, not u = (u_a+u_b)/2. So a 2:1 seam has
to be resolved by the TRIANGULATION - the hanging position becomes an ordinary
free node, and the coarse element is split so that it is a vertex on both sides.
The exported space is then a standard conforming P1 space on a graded mesh.

Checks, weakest to strongest:

  1. volume, boundary length, and int(n.x)dS == dim*volume  - geometry and the
     winding of the face elements.
  2. a Poisson solve runs and is sane.
  3. the LINEAR PATCH TEST: a linear field must be reproduced exactly. This is
     the standard check that a mesh has no spurious constraint or kink.
  4. CONFORMITY, which is the one that matters: every triangle edge must be
     shared by exactly two triangles or lie on the boundary. A quad split one way
     by one element and another way by its neighbour leaves dangling edges, and
     NO physics test can see it - a globally linear function lies in both
     triangulations, so the patch test passes at 1e-16 on a cracked mesh.
     toFinley() runs this internally (checkConformity) and throws if it fails,
     so simply getting a mesh back is the conformity check passing.

Run:
    ./bin/run-escript -n1 -t1 oxley-designs/validate_finley_export_hanging.py
"""
import sys

import numpy as np

import esys.escript as esc
import esys.finley as finley   # REQUIRED: registers the concrete domain type
import esys.oxley as oxley
from esys.escript import (ContinuousFunction, Function, FunctionOnBoundary,
                          Lsup, integrate, kronecker, whereZero)
from esys.escript.linearPDEs import LinearPDE, SolverOptions

rank = esc.getMPIRankWorld()
failures = []


def check(name, got, want, tol=1e-10):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    if rank == 0:
        print("    %-42s %18.12g (want %-16.12g) %s"
              % (name, got, want, "OK" if ok else "*** FAIL ***"))
    if not ok:
        failures.append(name)


def report(name, dom, lx, ly):
    """lx, ly: the domain extent, so the expected values are not hard-coded"""
    if rank == 0:
        print("\n%s" % name)
    fin = dom.toFinley()

    vol = integrate(esc.Scalar(1, Function(fin)))
    check("volume", vol, lx * ly)
    check("boundary length", integrate(esc.Scalar(1, FunctionOnBoundary(fin))),
          2 * (lx + ly))
    n = FunctionOnBoundary(fin).getNormal()
    xb = FunctionOnBoundary(fin).getX()
    check("int(n.x) dS   [face winding]", integrate(esc.inner(n, xb)), 2 * vol)

    # Per face, not just in the sum. int(n.x)dS is a global integral, so a
    # flipped face and a compensating one elsewhere would leave it intact. The
    # converter re-orients TRIANGLES by measuring their signed area, but passes
    # boundary faces through as the mesh view wound them, so nothing else would
    # catch a bad winding - and a face with an inward normal silently reverses
    # the Neumann/Robin term built on it.
    outward = esc.Vector(0, FunctionOnBoundary(fin))
    outward[0] = whereZero(xb[0] - lx) - whereZero(xb[0])
    outward[1] = whereZero(xb[1] - ly) - whereZero(xb[1])
    check("min(n . outward)  [every face]", esc.inf(esc.inner(n, outward)), 1.0)

    for tag, want in (("left", ly), ("right", ly), ("bottom", lx), ("top", lx)):
        mask = esc.Scalar(0, FunctionOnBoundary(fin))
        mask.setTaggedValue(tag, 1.0)
        check("boundary tag %-8s" % tag, integrate(mask), want)

    # linear patch test: exact reproduction of a linear field
    x = ContinuousFunction(fin).getX()
    u_exact = 1.0 + 2.0 * x[0] + 3.0 * x[1]
    pde = LinearPDE(fin, numEquations=1)
    pde.setValue(A=kronecker(fin))
    onbnd = (whereZero(x[0]) + whereZero(x[0] - lx)
             + whereZero(x[1]) + whereZero(x[1] - ly))
    pde.setValue(q=onbnd, r=u_exact)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
    pde.getSolverOptions().setTolerance(1e-12)
    u = pde.getSolution()
    check("linear patch test [rel Lsup]", Lsup(u - u_exact) / Lsup(u_exact), 0.0, 1e-9)

    # a Poisson solve must run and stay bounded
    pde2 = LinearPDE(fin, numEquations=1)
    pde2.setSymmetryOn()
    pde2.setValue(A=kronecker(fin), Y=1.0,
                  q=whereZero(x[0]) + whereZero(x[0] - lx))
    pde2.getSolverOptions().setTolerance(1e-12)
    # sup() is COLLECTIVE, so call it once on every rank and print the result
    # inside the guard. Calling it inside `if rank == 0` makes rank 0 issue one
    # more reduction than the others, and the run deadlocks at some later
    # collective, far from the cause.
    sup_s = esc.sup(pde2.getSolution())
    if rank == 0:
        print("    %-42s sup(u) = %.10g" % ("Poisson solve runs", sup_s))
    if not np.isfinite(sup_s):
        failures.append("Poisson solve")


# Cases, chosen from a survey of what balance actually produces. Blocks are
# always unit squares (l0=n0, l1=n1) so the elements are square.
#
# Note what the survey showed: a big level jump does NOT give a cell many hanging
# edges - balance turns it into a gradual grading, and 3 or 4 hanging edges arise
# only from an ISOLATED COARSE CELL ringed by finer ones. Both extremes are here.
CASES = [
    # name, per-block levels, (hanging edges per element seen in the survey)
    ("one seam [[2],[3]]",              [[2], [3]],                        "1"),
    ("cascade [[3,1],[1,2]]",           [[3, 1], [1, 2]],                  "1,2"),
    ("3x3 mixed",                       [[3,1,2],[1,2,1],[2,1,3]],         "1,2,3"),
    ("isolated coarse cell",            [[1,1,1],[1,0,1],[1,1,1]],         "4"),
    ("big jump [[0,3],[3,0]]",          [[0, 3], [3, 0]],                  "1,2"),
    ("diagonal ramp 4x4",               [[0,1,2,3],[1,2,3,2],[2,3,2,1],[3,2,1,0]], "1,2"),
    ("centre peak 5x5",                 [[0,0,1,0,0],[0,1,2,1,0],[1,2,3,2,1],
                                         [0,1,2,1,0],[0,0,1,0,0]],         "1,2"),
    ("checkerboard 4x4",                [[1,2,1,2],[2,1,2,1],[1,2,1,2],[2,1,2,1]], "1,2"),
    ("corner spike 3x3 (5 levels)",     [[4,0,0],[0,0,0],[0,0,2]],         "1,2"),
    ("stripes 4x3",                     [[3,0,3],[0,3,0],[3,0,3],[0,3,0]], "1,2"),
    ("control: uniform level 2",        2,                                 "none"),
]

for name, levels, hang in CASES:
    if isinstance(levels, int):
        n0 = n1 = 2
    else:
        n0, n1 = len(levels), len(levels[0])
    report("%s   [hanging edges: %s]" % (name, hang),
           oxley.Rectangle(n0=n0, n1=n1, l0=float(n0), l1=float(n1),
                           refine_level=levels),
           float(n0), float(n1))

if rank == 0:
    print("\nRESULT: %s" % ("ALL CHECKS PASSED" if not failures
                            else "FAILED (%s)" % ", ".join(failures)))
sys.exit(0 if not failures else 1)
