"""
A4 verification (3D): serial Poisson on a uniform conforming Block, solved on the
lnodes-based oxley assembly and compared against ripley on a matching grid.

Problem:  -div(grad u) = 1 on the unit cube, u = 0 on the boundary.
Compares order-independent scalars (sup u, integral of u).
"""
from esys.escript import (Solution, Function, sup, inf, integrate,
                          whereZero, kronecker)
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
import esys.oxley as oxley
import esys.ripley as ripley


def solve_poisson(dom):
    x = Solution(dom).getX()
    pde = LinearSinglePDE(dom)
    pde.setValue(A=kronecker(dom), Y=1.,
                 q=whereZero(x[0]) + whereZero(x[0] - 1.)
                 + whereZero(x[1]) + whereZero(x[1] - 1.)
                 + whereZero(x[2]) + whereZero(x[2] - 1.),
                 r=0.)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
    pde.getSolverOptions().setTolerance(1e-10)
    pde.getSolverOptions().setVerbosity(False)
    return pde.getSolution()


def report(name, u, dom):
    s, i = sup(u), inf(u)
    integ = integrate(u, Function(dom))
    print("  %-8s  sup=%.8f  inf=%.8f  integral=%.8f" % (name, s, i, integ))
    return s, i, integ


N = 8   # elements per axis (must match between domains)
print("3D Poisson  -Lap u = 1, u=0 on boundary,  %dx%dx%d elements" % (N, N, N))

# oxley: 2x2x2 blocks, refine_level=2  ->  2*2^2 = 8 elements per axis
odom = oxley.Block(numBlocks=(2, 2, 2), length=(1., 1., 1.), origin=(0., 0., 0.),
                   refine_level=2)
so = report("oxley", solve_poisson(odom), odom)

rdom = ripley.Brick(N, N, N, l0=1., l1=1., l2=1.)
sr = report("ripley", solve_poisson(rdom), rdom)

ds = abs(so[0] - sr[0])
dii = abs(so[1] - sr[1])
dint = abs(so[2] - sr[2])
print("  |dsup|=%.3e  |dinf|=%.3e  |dintegral|=%.3e" % (ds, dii, dint))
tol = 1e-6
print("RESULT:", "PASS" if (ds < tol and dii < tol and dint < tol) else "FAIL")
