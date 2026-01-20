
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
"""
example to demonstrate the use of Dirac Delta functions 
"""

__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"


from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE
from esys.weipa import saveVTK

try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
# generate domain:
if not HAVE_FINLEY:
    print("Finley module not available")
else:
    mydomain=Rectangle(30,30, l0=3, l1=2, 
                diracPoints=[(1.,1.), (2.,1.)],  diracTags=['in', 'out'])
    # fix the solution on the boundary
    x = mydomain.getX()
    gammaD = whereZero(x[0])+whereZero(x[1])+whereZero(x[0]-3.)+whereZero(x[1]-2.)
    # fix the solution on the boundary
    s=Scalar(0., DiracDeltaFunctions(mydomain))
    s.setTaggedValue('in', +1.)
    s.setTaggedValue('out', -1.)
    # define PDE and get its solution u
    mypde = LinearSinglePDE(domain=mydomain)
    mypde.setValue(q=gammaD, A=kronecker(2), y_dirac=s)
    u = mypde.getSolution()
    print("Solution = ",str(u))
    # write u to an external file
    saveVTK("output/u.vtu",sol=u)
