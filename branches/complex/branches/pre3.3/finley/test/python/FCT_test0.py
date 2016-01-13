
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import TransportPDE
from esys.finley import Rectangle, Brick
#from esys.ripley import Rectangle, Brick
from esys.weipa import saveVTK
from math import pi, ceil

NE=50
dom=Rectangle(NE,1,l1=1./NE)
dom=Rectangle(NE,NE)
fc=TransportPDE(dom,numEquations=1)
fc.getSolverOptions().setVerbosityOn()
fc.getSolverOptions().setODESolver(fc.getSolverOptions().LINEAR_CRANK_NICOLSON)
fc.getSolverOptions().setODESolver(fc.getSolverOptions().BACKWARD_EULER)
fc.getSolverOptions().setODESolver(fc.getSolverOptions().CRANK_NICOLSON)
fc.setValue(M=1,C=[-1,0])
x=dom.getX()
u0=whereNegative(x[0]-1./NE)

c=0
t=0

saveVTK("u.%s.vtu"%c,u=u0)
fc.setInitialSolution(u0)
dt=fc.getSafeTimeStepSize() 

print "u0 =",u0
T_END=dt
print "dt = ",dt
while t<T_END:
    print("time step t=",t+dt)	
    u=fc.getSolution(dt)
    saveVTK("u.%s.vtu"%(c+1,),u=u)
    print "u =",u
    c+=1
    t+=dt
