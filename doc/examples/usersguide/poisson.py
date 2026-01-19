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
__url__="https://github.com/LutzGross/esys-escript.github.io"

from esys.escript import *
from esys.escript.linearPDEs import Poisson
try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
from esys.weipa import saveVTK
# generate domain:
if not HAVE_FINLEY:
    print("Finley module not available")
else:
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

