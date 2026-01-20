##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

from esys.escript import *
from esys.escript.linearPDEs import Poisson
try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
from esys.weipa import saveVTK

if not HAVE_FINLEY:
    print("Finley module not available")
else:
    # generate domain:
    mydomain = Rectangle(l0=1.,l1=1.,n0=40, n1=20)
    # define characteristic function of Gamma^D
    x = mydomain.getX()
    gammaD = whereZero(x[0])+whereZero(x[1])
    # define PDE and get its solution u
    mypde = Poisson(domain=mydomain)
    mypde.setValue(f=1,q=gammaD)
    u = mypde.getSolution()
    # write u to an external file
    saveVTK("output/u.vtu",sol=u)

