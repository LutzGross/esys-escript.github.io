# 3D point-source Poisson: dirac point in the domain; check the solve responds
from esys.escript import Solution, DiracDeltaFunctions, Scalar, whereZero, kronecker, Lsup, sup, integrate, Function
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
import esys.oxley as oxley, esys.ripley as ripley
def run(name, mk):
    d = mk([(0.5,0.5,0.5)], ["src"])
    x = Solution(d).getX()
    pde = LinearSinglePDE(d)
    y_dirac = Scalar(0., DiracDeltaFunctions(d)); y_dirac.setTaggedValue("src", 1.0)
    q = whereZero(x[0])+whereZero(x[0]-1)+whereZero(x[1])+whereZero(x[1]-1)+whereZero(x[2])+whereZero(x[2]-1)
    pde.setValue(A=kronecker(d), y_dirac=y_dirac, q=q)
    pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
    u = pde.getSolution()
    print("  %-8s sup(u)=%.6f  int(u)=%.6f"%(name, sup(u), integrate(u,Function(d))))
run("oxley",  lambda p,t: oxley.Brick(n0=8,n1=8,n2=8, diracPoints=p, diracTags=t))
run("ripley", lambda p,t: ripley.Brick(8,8,8, diracPoints=p, diracTags=t))
