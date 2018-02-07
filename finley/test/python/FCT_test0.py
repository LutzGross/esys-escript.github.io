
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"


from esys.escript import *
from esys.escript.linearPDEs import TransportPDE, SolverOptions
from esys.finley import Rectangle, Brick
#from esys.ripley import Rectangle, Brick
from esys.weipa import saveVTK
from math import pi, ceil

NE=50
dom=Rectangle(NE,1,l1=1./NE)
dom=Rectangle(NE,NE)
fc=TransportPDE(dom,numEquations=1)
fc.getSolverOptions().setVerbosityOn()
fc.getSolverOptions().setODESolver(SolverOptions.LINEAR_CRANK_NICOLSON)
fc.getSolverOptions().setODESolver(SolverOptions.BACKWARD_EULER)
fc.getSolverOptions().setODESolver(SolverOptions.CRANK_NICOLSON)
fc.setValue(M=1,C=[-1,0])
x=dom.getX()
u0=whereNegative(x[0]-1./NE)

c=0
t=0

saveVTK("u.%s.vtu"%c,u=u0)
fc.setInitialSolution(u0)
dt=fc.getSafeTimeStepSize() 

print("u0 =",u0)
T_END=dt
print("dt = ",dt)
while t<T_END:
    print("time step t=",t+dt)
    u=fc.getSolution(dt)
    saveVTK("u.%s.vtu"%(c+1,),u=u)
    print("u =",u)
    c+=1
    t+=dt
