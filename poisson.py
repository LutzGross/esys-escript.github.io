from esys.escript import *
from esys.escript.linearPDEs import LinearPDE,SolverOptions
from esys.ripley import Rectangle, Brick
from time import time

BLOCKSIZE=1

dom = Rectangle(l0=1.,l1=1.,n0=599, n1=599)
#dom = Brick(l0=1.,l1=1.,n0=51, n1=51, n2=17)
x = dom.getX()
gammaD = whereZero(x[0])+whereZero(x[1])

def solve(target):
    pde = LinearPDE(dom, numEquations=BLOCKSIZE, numSolutions=BLOCKSIZE)
    A = pde.createCoefficient("A")
    q = pde.createCoefficient("q")
    Y_reduced = pde.createCoefficient("Y_reduced")

    if BLOCKSIZE == 1:
        A = kronecker(dom)
        q = gammaD
        Y_reduced = 1.
    else:
        for i in range(BLOCKSIZE):
            A[i,:,i,:] = kronecker(dom)
            q[i] = gammaD
            Y_reduced[i] = 1.

    pde.setValue(A=A, Y_reduced=Y_reduced, q=q)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
    pde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
    pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
    pde.getSolverOptions().setSolverTarget(target)
    pde.getSolverOptions().setPreconditioner(SolverOptions.NO_PRECONDITIONER)
    pde.getSolverOptions().setVerbosityOn()
    #pde.setDebugOn()
    #rhs=pde.getRightHandSide()
    #saveDataCSV('/tmp/rhs.csv',rhs=rhs)
    #pde.getSystem()[0].saveMM('/tmp/poissonripley.mtx')
    t0=time()
    print("pde.getSolution()")
    x = pde.getSolution()
    t1=time()
    print("Solver Time: %s"%(t1-t0))
    print("Solution: %s..%s"%(inf(x),sup(x)))
    print(x)


for target in SolverOptions.TARGET_CPU,SolverOptions.TARGET_GPU:
    solve(target)
    print("-------")

