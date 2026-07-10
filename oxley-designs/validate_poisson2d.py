"""
A4 verification: serial 2D Poisson on a uniform conforming Block, solved on the
lnodes-based oxley assembly and compared against ripley on a matching grid.

Problem:  -div(grad u) = 1 on the unit square, u = 0 on the boundary.
We compare order-independent scalars (sup u, integral of u) so node/DOF
orderings need not match between the two domains.
"""
import numpy as np
from esys.escript import (Solution, Function, ContinuousFunction,
                          sup, inf, integrate, whereZero, Lsup)
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator
import esys.oxley as oxley
import esys.ripley as ripley


def solve_poisson(dom):
    x = Solution(dom).getX()
    pde = LinearSinglePDE(dom)
    pde.setValue(A=kronecker_like(dom),
                 Y=1.,
                 q=whereZero(x[0]) + whereZero(x[0] - 1.)
                 + whereZero(x[1]) + whereZero(x[1] - 1.),
                 r=0.)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
    pde.getSolverOptions().setTolerance(1e-10)
    pde.getSolverOptions().setVerbosity(False)
    return pde.getSolution()


def kronecker_like(dom):
    from esys.escript import kronecker
    return kronecker(dom)


def report(name, u, dom):
    s = sup(u)
    i = inf(u)
    integ = integrate(u, Function(dom))
    print("  %-8s  sup=%.8f  inf=%.8f  integral=%.8f" % (name, s, i, integ))
    return s, i, integ


N = 16   # elements per axis (must match between domains)

print("2D Poisson  -Lap u = 1, u=0 on boundary,  %dx%d elements" % (N, N))

# oxley: 4x4 blocks, refine_level=2  ->  4*2^2 = 16 elements per axis
odom = oxley.Block(numBlocks=(4, 4), length=(1., 1.), origin=(0., 0.),
                   refine_level=2)
uo = solve_poisson(odom)
so = report("oxley", uo, odom)

# ripley: 16x16 element grid over the unit square
rdom = ripley.Rectangle(N, N, l0=1., l1=1.)
ur = solve_poisson(rdom)
sr = report("ripley", ur, rdom)

ds = abs(so[0] - sr[0])
di = abs(so[2] - sr[2])
print("  |dsup|=%.3e  |dintegral|=%.3e" % (ds, di))

tol = 1e-6
ok = ds < tol and di < tol
print("RESULT:", "PASS" if ok else "FAIL")
