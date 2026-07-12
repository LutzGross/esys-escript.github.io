# 3D linear elasticity: -div(C:grad u) = 0, Dirichlet u; A couples components.
from esys.escript import Solution, sup, inf, whereZero, kronecker, Lsup, Function, integrate
from esys.escript.linearPDEs import LinearPDE, SolverOptions
import esys.oxley as oxley, esys.ripley as ripley
def solve(dom):
    DIM=3
    lam, mu = 1.0, 1.0
    pde = LinearPDE(dom, numEquations=DIM, numSolutions=DIM)
    A = pde.createCoefficient("A")   # A[i,j,k,l] = lam d_ij d_kl + mu(d_ik d_jl + d_il d_jk)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    v=0.
                    if i==j and k==l: v+=lam
                    if i==k and j==l: v+=mu
                    if i==l and j==k: v+=mu
                    A[i,j,k,l]=v
    x = Solution(dom).getX()
    # Dirichlet: u = (x0, 0, 0)*0.1 on all boundaries (rigid stretch), interior solves
    onb = whereZero(x[0])+whereZero(x[0]-1)+whereZero(x[1])+whereZero(x[1]-1)+whereZero(x[2])+whereZero(x[2]-1)
    q = pde.createCoefficient("q"); r = pde.createCoefficient("r")
    for i in range(DIM): q[i]=onb
    r[0]=0.1*x[0]
    pde.setValue(A=A, q=q, r=r)
    pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
    return pde.getSolution()
for n,dom in [("ripley3D",ripley.Brick(8,8,8)),("oxley3D",oxley.Block(numBlocks=(2,2,2),refine_level=2))]:
    try:
        u=solve(dom)
        print("  %-9s u0[%.6f,%.6f] u1[%.6f,%.6f] u2[%.6f,%.6f] |u|int=%.6f"%(n,
            inf(u[0]),sup(u[0]),inf(u[1]),sup(u[1]),inf(u[2]),sup(u[2]), integrate(Lsup(u)*0+u[0],Function(dom))))
    except Exception as ex:
        print("  %-9s FAIL %s"%(n,repr(ex)[:90]))
