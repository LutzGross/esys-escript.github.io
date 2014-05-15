from esys.escript import *
from esys.escript.linearPDEs import LinearPDE,SolverOptions
from esys.ripley import Rectangle
from time import time

BLOCKSIZE=4

dom = Rectangle(l0=1.,l1=1.,n0=2, n1=2)
x = dom.getX()
n = dom.getNormal()

pde = LinearPDE(dom, numEquations=BLOCKSIZE, numSolutions=BLOCKSIZE)
A = pde.createCoefficient("A")
D = pde.createCoefficient("D")
Y = pde.createCoefficient("Y")
d = pde.createCoefficient("d")
y = pde.createCoefficient("y")
if BLOCKSIZE == 1:
    omega=0.1
    A[:,:]=kronecker(dom)
    D=omega
    Y=omega*x[0]
    d=10
    y=n[0]+10*x[0]
else:
    for i in range(BLOCKSIZE):
        omega=0.1*i
        A[i,:,i,:]=kronecker(dom)
        D[i,:]=[omega]*BLOCKSIZE
        Y[i]=omega*x[0]
        d[i,:]=[10]*BLOCKSIZE
        y[i]=n[0]+10*x[0]

pde.setValue(A=A, D=D, Y=Y, d=d, y=y)
pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
pde.getSolverOptions().setPreconditioner(SolverOptions.NO_PRECONDITIONER)
pde.getSolverOptions().setVerbosityOn()
#pde.setDebugOn()
rhs=pde.getRightHandSide()
#saveDataCSV('helmholtz_rhs.csv',rhs=rhs)
pde.getSystem()[0].saveMM('/tmp/helmholtzripley.mtx')
t0=time()
x = pde.getSolution()
t1=time()
print "Paso Solver Time: ", t1-t0
print "Solution: %s..%s"%(inf(x),sup(x))
print x

