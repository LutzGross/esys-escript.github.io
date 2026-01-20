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

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://github.com/LutzGross/esys-escript.github.io"

from esys.escript import *
from esys.escript.pdetools import MaskFromBoundaryTag
from esys.escript.linearPDEs import LinearPDE

try:
    from esys.finley import ReadGmsh
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
from esys.weipa import saveVTK

if not HAVE_FINLEY:
    print("Finley module not available")
else:
    # ... generate domain ...
    print("read in mesh")
    domain = ReadGmsh("inclusion.msh", 2, optimize=True)
    # ... set some parameters ...
    kappa = 1.
    q_c = 0.1
    q_out = 0.05
    T_ref = 15.

    # set source term:
    q_H=Scalar(0., Function(domain))
    q_H.setTaggedValue("Inclusion", q_c)
    # set flux term:
    q_S=Scalar(0., FunctionOnBoundary(domain))
    q_S.setTaggedValue("Left", q_out)
    q_S.setTaggedValue("Right", q_out)
    # ... set temperature at the bottom of the domain:
    # y=domain.getX()[1]
    # fixedT1 = whereZero(y - inf(y))
    fixedT = MaskFromBoundaryTag(domain, *('Bottom', ))
    # ... open PDE and set coefficients ...
    mypde = LinearPDE(domain)
    mypde.setSymmetryOn()
    mypde.setValue(A=kappa * kronecker(domain), Y=q_H, y=q_S, \
                   q=fixedT, r=T_ref)
    # ... calculate error of the PDE solution ...
    T = mypde.getSolution()
    print("Temperature is ", T)
    saveVTK("output/T.vtu", temperature=T)

