from esys.escript import *
from esys.escript.linearPDEs import Poisson,SolverOptions
from esys.ripley import Rectangle
from time import time

mydomain = Rectangle(l0=1.,l1=1.,n0=9, n1=9)
x = mydomain.getX()
gammaD = whereZero(x[0])+whereZero(x[1])
mypde = Poisson(domain=mydomain)
#mypde.getSolverOptions().setPackage(SolverOptions.CUSP)
mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
mypde.getSolverOptions().setPreconditioner(SolverOptions.NO_PRECONDITIONER)
mypde.getSolverOptions().setVerbosityOn()
#mypde.setDebugOn()
mypde.setValue(f_reduced=1,q=gammaD)
#rhs=mypde.getRightHandSide()
#saveDataCSV('/tmp/rhs.csv',rhs=rhs)
#mypde.getSystem()[0].saveMM('/tmp/poisson.mtx')
#print rhs

t0=time()
x = mypde.getSolution()
t1=time()
print "Time: ", t1-t0
print x

