##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
from esys.weipa import saveVTK

if not HAVE_FINLEY:
    print("Finley module not available")
else:
    #... set some parameters ...
    kappa=1.
    omega=0.1
    eta=10.
    #... generate domain ...
    mydomain = Rectangle(l0=5.,l1=1.,n0=50, n1=10)
    #... open PDE and set coefficients ...
    mypde=LinearPDE(mydomain)
    mypde.setSymmetryOn()
    n=mydomain.getNormal()
    x=mydomain.getX()
    mypde.setValue(A=kappa*kronecker(mydomain),D=omega,Y=omega*x[0], \
                   d=eta,y=kappa*n[0]+eta*x[0])
    #... calculate error of the PDE solution ...
    u=mypde.getSolution()
    print("error is ",Lsup(u-x[0]))
    # output should be similar to "error is 1.e-7"
    saveVTK("x0.vtu",sol=u)
 
