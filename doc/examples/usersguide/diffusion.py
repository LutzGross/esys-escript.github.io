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
#... set some parameters ...
if not HAVE_FINLEY:
    print("Finley module not available")
else:
    xc=[0.02,0.002]
    r=0.001
    qc=50.e6
    Tref=0.
    rhocp=2.6e6
    eta=75.
    kappa=240.
    tend=5.
    # ... time, time step size and counter ...
    t=0
    h=0.1
    i=0
    #... generate domain ...
    mydomain = Rectangle(l0=0.05,l1=0.01,n0=250, n1=50)
    #... open PDE ...
    mypde=LinearPDE(mydomain)
    mypde.setSymmetryOn()
    mypde.setValue(A=kappa*kronecker(mydomain),D=rhocp/h,d=eta,y=eta*Tref)
    # ... set heat source: ....
    x=mydomain.getX()
    qH=qc*whereNegative(length(x-xc)-r)
    # ... set initial temperature ....
    T=Tref
    # ... start iteration:
    while t<tend:
          i+=1
          t+=h
          print("time step :",t)
          mypde.setValue(Y=qH+rhocp/h*T)
          T=mypde.getSolution()
          saveVTK("output/T.%d.vtu"%i,temp=T)

