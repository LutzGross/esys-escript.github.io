from __future__ import division
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
    # set dimensions
    L0=1.;L1=1. 
    # bottom temperature:
    T_bot=100
    # location, size and value of heat source
    xc=[0.3,0.4]; r=0.1; Qc=3000
    # create domain
    mydomain=Rectangle(l0=L0,l1=L1,n0=20,n1=20)
    x=mydomain.getX()
    k=1
    # temperature for boundary condition
    T_D=T_bot/L1*(L1-x[1])
    # heat source
    Q=Qc*whereNegative(length(x-xc)-r)
    # create PDE:
    mypde=LinearPDE(mydomain)
    mypde.setSymmetryOn()
    # set coefficients:
    mypde.setValue(A=k*kronecker(mydomain),Y=Q, r=T_D, \
                    q=whereZero(x[1])+whereZero(x[1]-L1))
    # get temperature:
    T=mypde.getSolution()
    # write to file:
    saveVTK("u.vtu",T=T)

