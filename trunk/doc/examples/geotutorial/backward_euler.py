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
from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

# import tools
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
try:
    from esys.dudley import Rectangle
    HAVE_DUDLEY = True
except ImportError:
    HAVE_DUDLEY = False
from esys.weipa import saveVTK

if not HAVE_DUDLEY:
    print("Dudley module not available")
else:
    # end of simulation time
    t_end=0.1
    # time step size:
    dt=0.01
    # dimensions:
    L0=1.;L1=1.
    # location, size and value of heat source
    xc=[0.3,0.4]; r=0.1; Qc=3000
    # material parameter
    k=1; rhocp=100; 
    # bottom temperature:
    T_bot=100
    # generate domain:
    mydomain=Rectangle(l0=L0,l1=L1,n0=20,n1=20)
    x=mydomain.getX()
    # set boundray temperature:
    T_D=T_bot/L1*(L1-x[1])
    # set heat source:
    Q=Qc*whereNegative(length(x-xc)-r)
    # generate domain:
    mypde=LinearPDE(mydomain)
    mypde.setSymmetryOn()
    # set PDE coefficients:
    mypde.setValue(A=dt*k*kronecker(mydomain), D=dt*rhocp, 
                    r=T_D, q=whereZero(x[1])+whereZero(x[1]-L1))
    # initial temperature
    T=T_D 
    # step counter and time marker:
    N=0; t=0
    # stop when t_end is reached:
    while t<t_end:
        print("time step %d, t=%s"%(N,t))
        # update PDE coefficient:
        mypde.setValue(Y=dt*rhocp*T+dt*Q)
        # new temperature:
        T=mypde.getSolution()
        # save as VTK for visualisation:
        saveVTK("u.%s.vtu"%N,T=T)
        # increase counter and marker:
        N+=1; t+=dt
