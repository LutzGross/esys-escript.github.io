from esys.escript import Solution, FunctionOnBoundary, sup, inf, kronecker, integrate, Function
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
import esys.oxley as oxley, esys.ripley as ripley
def solve(d):
    xb = FunctionOnBoundary(d).getX()
    pde = LinearSinglePDE(d)
    pde.setValue(A=kronecker(d), d=1., y=xb[0])   # -Lap u=0, Robin: n.grad u + u = x0
    pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
    return pde.getSolution()
for n,dom in [("oxley3D",oxley.Block(numBlocks=(2,2,2),refine_level=2)),("ripley3D",ripley.Brick(8,8,8))]:
    u=solve(dom); print("  %-9s sup=%.7f inf=%.7f  int_u=%.7f"%(n,sup(u),inf(u),integrate(u,Function(dom))))
# also test constant Neumann-ish: d=1, y=2 -> u=2 everywhere (interior Lap=0, bc u=2)
print("--- d=1,y=2 (expect u=2) ---")
for n,dom in [("oxley3D",oxley.Block(numBlocks=(2,2,2),refine_level=2)),("ripley3D",ripley.Brick(8,8,8))]:
    pde=LinearSinglePDE(dom); pde.setValue(A=kronecker(dom), d=1., y=2.)
    pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
    u=pde.getSolution(); print("  %-9s sup=%.7f inf=%.7f"%(n,sup(u),inf(u)))
