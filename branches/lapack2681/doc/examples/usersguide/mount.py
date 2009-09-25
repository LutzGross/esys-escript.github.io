########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"
from esys.escript import *
from esys.escript.models import Mountains
from esys.finley import Brick,Rectangle
from math import pi, ceil

NE=16
DIM=3
H=1.
L=2*H
OMEGA=10
EPS=0.01
t=0
T_END=0.05 # set T_END=(2*pi)/OMEGA to run a full simulation
n=0
if DIM==2:
  mydomain=Rectangle(int(ceil(L*NE/H)),NE,l0=L,l1=H,order=1, useFullElementOrder=True,optimize=True)
else:
  mydomain=Brick(int(ceil(L*NE/H)),int(ceil(L*NE/H)),NE,l0=L,l1=L,l2=H,order=1, useFullElementOrder=True,optimize=True)

x=mydomain.getX()
v = Vector(0.0, Solution(mydomain))
if DIM==2:
    a0=1
    n0=1
    n1=0.5
    a1=-(a0*n0)/n1
    v[0]=a0*sin(pi*n0*x[0])* cos(pi*n1*x[1])
    v[1]=a1*cos(pi*n0*x[0])* sin(pi*n1*x[1])
else:
    a0=1
    a1=1
    n0=2
    n1=2
    n2=0.5
    a2=-(a0*n0+a1*n1)/n2
    v[0]=a0*sin(pi*n0*x[0])* cos(pi*n1*x[1])* cos(pi*n2*x[2])
    v[1]=a1*cos(pi*n0*x[0])* sin(pi*n1*x[1])* cos(pi*n2*x[2])
    v[2]=a2*cos(pi*n0*x[0])* cos(pi*n1*x[1])* sin(pi*n2*x[2])


mts=Mountains(mydomain,eps=EPS)
while t<T_END:
    print "STEP ", t
    mts.setVelocity(v*cos(OMEGA*t))
    Z=mts.update()
    
    saveVTK("state.%d.vtu"%n,sol=Z, v=mts.getVelocity())
    print "Integral(Z)=",integrate(Z),Lsup(mts.getVelocity()[DIM-1])
    n+=1
    t+=mts.getSafeTimeStepSize()

