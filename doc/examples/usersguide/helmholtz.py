
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

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
    saveVTK("output/x0.vtu",sol=u)
 
