from __future__ import division, print_function
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

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"
from esys.escript import *
from esys.escript.models import Mountains
from esys.weipa import saveVTK
from math import pi, ceil

try:
    from esys.finley import Rectangle, Brick
    HAVE_FINLEY = True
except ImportError:
    print("Finley module required but not available")
    HAVE_FINLEY = False

if HAVE_FINLEY:
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
        print("STEP ", t)
        mts.setVelocity(v*cos(OMEGA*t))
        Z=mts.update()
        
        saveVTK("state.%d.vtu"%n,sol=Z, v=mts.getVelocity())
        print("Integral(Z)=",integrate(Z),Lsup(mts.getVelocity()[DIM-1]))
        n+=1
        t+=mts.getSafeTimeStepSize()

