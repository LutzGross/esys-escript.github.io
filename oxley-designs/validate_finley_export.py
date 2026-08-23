"""
Spike validation: oxley forest -> finley mesh, conforming and graded.

Checks, in increasing order of how much of the pipeline they exercise:
  1. volume           - integrate(1) over Function
  2. boundary area    - integrate(1) over FunctionOnBoundary
  3. NORMALS          - integrate(n.x) dS == dim * volume. This is the one that
                        catches a face element wound the wrong way round, which
                        is otherwise silent and corrupts Neumann/Robin BCs.
  4. tags             - each boundary tag has the area it should
  5. patch test       - a linear field solved exactly
  6. Poisson solve    - against the equivalent ripley domain

The 3D cases run twice: on a conforming forest, and on the GRADED forests from
the shared table, where a coarse octant is cut into up to forty-eight
tetrahedra and a boundary quad into up to six triangles.

Run FROM THE PROJECT ROOT, which is where the shared table of forests is found:
    ./bin/run-escript -n1 oxley-designs/validate_finley_export.py
    ./bin/run-escript -n3 oxley-designs/validate_finley_export.py

Under MPI compare oxley against oxley SERIAL, not against ripley: ripley adjusts
its own mesh to decompose across ranks ("Domain setup has been adjusted"), so the
ripley reference moves with the rank count while oxley's does not. At 4 ranks
ripley gives 0.123456790123 where oxley->finley still gives exactly 0.125.

TRAP, cost an hour here: integrate(), sup(), inf(), Lsup() and isConforming() are
all COLLECTIVE. Never call one inside an `if rank == 0:` block - rank 0's
reduction then pairs with whatever collective the other ranks reached next, the
ranks drift out of step and the run deadlocks somewhere later, with the failure
nowhere near the cause. Evaluate on every rank, then print on rank 0.
"""
import sys

sys.path.insert(0, "oxley/test/python")   # for the shared table of forests

from esys.escript import (Function, FunctionOnBoundary, Solution,
                          ContinuousFunction, Lsup, integrate, inf, sup,
                          whereZero, grad)
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.util import interpolate
import esys.finley as finley
import esys.oxley as oxley
import esys.ripley as ripley
import esys.escript as esc
import oxley_meshes

TOL = 1e-10
failures = []
rank = esc.getMPIRankWorld()
nranks = esc.getMPISizeWorld()


def check(name, got, want, tol=TOL):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    if rank == 0:
        print("  %-38s %20.12g  (want %-18.12g) %s"
              % (name, got, want, "OK" if ok else "*** FAIL ***"))
    if not ok:
        failures.append(name)
    return ok


def check_mesh(dom, dim, volume, area):
    """geometry and, crucially, normal orientation"""
    check("integrate(1) over Function", integrate(esc.Scalar(1, Function(dom))),
          volume)
    check("integrate(1) over boundary",
          integrate(esc.Scalar(1, FunctionOnBoundary(dom))), area)

    # divergence theorem: int_dOmega n.x dS = int_Omega div(x) dV = dim*volume.
    # A face wound the wrong way flips its normal and this drops immediately.
    fob = FunctionOnBoundary(dom)
    n = dom.getNormal()
    x = fob.getX()
    check("integrate(n.x) dS  [normals]", integrate(esc.inner(n, x)),
          dim * volume)


def poisson(dom, tag_low, tag_high):
    """-div(grad(u)) = 1 with u=0 on two opposite faces"""
    pde = LinearPDE(dom, numEquations=1)
    pde.setSymmetryOn()
    x = dom.getX()
    pde.setValue(A=esc.kronecker(dom), Y=1.0,
                 q=whereZero(x[0]) + whereZero(x[0] - sup(x[0])))
    pde.getSolverOptions().setTolerance(1e-12)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
    return pde.getSolution()


def patch_test(dom, dim):
    """Solve with a linear exact solution imposed all around the boundary.

    Linear fields are in the Tet3/Tri3 space exactly, so any error here is a
    defect in the mesh: a flipped element, a bad coordinate, a wrong node
    reference. NOTE it does NOT detect a quad face split one way by one element
    and the other way by its neighbour - a globally linear function lies in both
    triangulations of a planar quad. That crack is caught combinatorially, by
    the face-hash check inside toFinley().
    """
    x = dom.getX()
    u_ex = 1.0 + 2.0 * x[0] + 3.0 * x[1]
    if dim == 3:
        u_ex = u_ex + 4.0 * x[2]
    on_bnd = (whereZero(x[0]) + whereZero(x[0] - 1.0)
              + whereZero(x[1]) + whereZero(x[1] - 1.0))
    if dim == 3:
        on_bnd = on_bnd + whereZero(x[2]) + whereZero(x[2] - 1.0)
    pde = LinearPDE(dom, numEquations=1)
    pde.setSymmetryOn()
    pde.setValue(A=esc.kronecker(dom), Y=0.0, q=on_bnd, r=u_ex)
    pde.getSolverOptions().setTolerance(1e-14)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
    u = pde.getSolution()
    return Lsup(u - u_ex) / Lsup(u_ex)


def run_2d():
    if rank == 0:
        print("\n=== 2D: oxley Rectangle(2,2, refine_level=2) -> finley ===")
    # 2 x 2 blocks, each refined twice -> 8 x 8 = 64 elements over [0,1]^2
    ox = oxley.Rectangle(n0=2, n1=2, l0=1.0, l1=1.0, refine_level=2)
    conforming = ox.isConforming()   # collective: every rank must call it
    if rank == 0:
        print("  forest is conforming:", conforming)
    fin = ox.toFinley()
    if rank == 0:
        print("  -> finley domain, dim =", fin.getDim())

    check_mesh(fin, 2, 1.0, 4.0)

    # each side of the unit square has length 1
    for tag, want in (("left", 1.0), ("right", 1.0),
                      ("bottom", 1.0), ("top", 1.0)):
        m = esc.Scalar(0, FunctionOnBoundary(fin))
        m.setTaggedValue(tag, 1.0)
        check("boundary tag '%s' length" % tag, integrate(m), want)

    check("linear patch test  [relative Lsup]", patch_test(fin, 2), 0.0, tol=1e-11)

    u = poisson(fin, "left", "right")
    rip = ripley.Rectangle(n0=8, n1=8, l0=1.0, l1=1.0)
    ur = poisson(rip, "left", "right")
    # evaluate the collectives on EVERY rank before printing on rank 0 only
    su, sur = sup(u), sup(ur)
    if rank == 0:
        print("  oxley->finley sup(u) = %.12g" % su)
        print("  ripley Rec4   sup(u) = %.12g" % sur)
    # exact solution of -u'' = 1, u(0)=u(1)=0 is x(1-x)/2, peak 1/8. Tri3 on this
    # mesh is only approximate, so this is a sanity band, not an identity.
    check("Poisson sup(u) near analytic 1/8", su, 0.125, tol=5e-2)


def run_3d():
    if rank == 0:
        print("\n=== 3D: oxley Brick(2,2,2, refine_level=1) -> finley ===")
    # 2^3 blocks each refined once -> 4 x 4 x 4 = 64 elements over [0,1]^3
    ox = oxley.Brick(n0=2, n1=2, n2=2, l0=1.0, l1=1.0, l2=1.0, refine_level=1)
    conforming = ox.isConforming()   # collective: every rank must call it
    if rank == 0:
        print("  forest is conforming:", conforming)
    fin = ox.toFinley()
    if rank == 0:
        print("  -> finley domain, dim =", fin.getDim())

    check_mesh(fin, 3, 1.0, 6.0)

    for tag, want in (("left", 1.0), ("right", 1.0), ("bottom", 1.0),
                      ("top", 1.0), ("front", 1.0), ("back", 1.0)):
        m = esc.Scalar(0, FunctionOnBoundary(fin))
        m.setTaggedValue(tag, 1.0)
        check("boundary tag '%s' area" % tag, integrate(m), want)

    check("linear patch test  [relative Lsup]", patch_test(fin, 3), 0.0, tol=1e-11)

    u = poisson(fin, "left", "right")
    su = sup(u)   # collective: every rank must call it
    if rank == 0:
        print("  oxley->finley sup(u) = %.12g" % su)
    check("Poisson sup(u) near analytic 1/8", su, 0.125, tol=5e-2)


def run_3d_graded():
    """
    The same checks on forests with 2:1 seams, which is where the 3D split has
    something to do.

    A coarse octant then meets finer ones, so its faces carry a centre node and
    its edges midpoints, it is cut into up to forty-eight tetrahedra coned from
    a node at its own centre, and a boundary quad becomes up to six triangles
    rather than two. The boundary checks are the ones that gain most from that:
    a triangle dropped from a polygon shows up in the area, and one wound the
    wrong way in the divergence theorem.

    What is NOT checked here, and cannot be: whether neighbouring octants cut
    their shared face the same way. A globally linear function lies in both
    triangulations of a planar quad, so the patch test below passes at machine
    precision on a mesh full of cracks, and so do the volume and the area. That
    one is caught combinatorially by the face hash inside toFinley(), which runs
    on every export - so reaching these checks at all already means it passed.
    """
    for name, levels in (("seam", oxley_meshes.SEAM_3D),
                         ("mixed", oxley_meshes.MIXED_3D),
                         ("isolated", oxley_meshes.ISOLATED_3D)):
        if rank == 0:
            print("\n=== 3D graded: oxley Brick, %s -> finley ===" % name)
        ox = oxley_meshes.graded3D(levels)
        conforming = ox.isConforming()   # collective: every rank must call it
        if rank == 0:
            print("  forest is conforming:", conforming)
        if conforming:
            # the case would silently become another conforming one
            failures.append("%s: the forest has no 2:1 seam" % name)
            if rank == 0:
                print("  *** FAIL *** this case is supposed to be graded")
            continue
        fin = ox.toFinley()

        check_mesh(fin, 3, 1.0, 6.0)

        for tag in ("left", "right", "bottom", "top", "front", "back"):
            m = esc.Scalar(0, FunctionOnBoundary(fin))
            m.setTaggedValue(tag, 1.0)
            check("boundary tag '%s' area" % tag, integrate(m), 1.0)

        check("linear patch test  [relative Lsup]", patch_test(fin, 3), 0.0,
              tol=1e-11)

        u = poisson(fin, "left", "right")
        su = sup(u)   # collective: every rank must call it
        if rank == 0:
            print("  oxley->finley sup(u) = %.12g" % su)
        check("Poisson sup(u) near analytic 1/8", su, 0.125, tol=5e-2)


if __name__ == "__main__":
    if rank == 0:
        print("ranks:", nranks)
    run_2d()
    run_3d()
    run_3d_graded()
    if rank == 0:
        print()
        if failures:
            print("FAILED %d check(s): %s" % (len(failures), ", ".join(failures)))
        else:
            print("ALL CHECKS PASSED")
    if failures:
        sys.exit(1)
