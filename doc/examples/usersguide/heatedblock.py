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
from esys.weipa import saveVTK

try:
    from esys.finley import Brick
    HAVE_FINLEY = True
except ImportError:
    print("Finley module required but not available")
    HAVE_FINLEY = False

if HAVE_FINLEY:

    #... set some parameters ...
    lam=1.
    mu=0.1
    alpha=1.e-6
    xc=[0.3,0.3,1.]
    beta=8.
    T_ref=0.
    T_0=1.
    #... generate domain ...
    mydomain = Brick(l0=1.,l1=1., l2=1.,n0=10, n1=10, n2=10)
    x=mydomain.getX()
    #... set temperature ...
    T=T_0*exp(-beta*length(x-xc))
    #... open symmetric PDE ...
    mypde=LinearPDE(mydomain)
    mypde.setSymmetryOn()
    #... set coefficients ...
    C=Tensor4(0.,Function(mydomain))
    for i in range(mydomain.getDim()):
      for j in range(mydomain.getDim()):
         C[i,i,j,j]+=lam
         C[i,j,i,j]+=mu
         C[i,j,j,i]+=mu
    msk=whereZero(x[0])*[1.,0.,0.] \
       +whereZero(x[1])*[0.,1.,0.] \
       +whereZero(x[2])*[0.,0.,1.]
    sigma0=(lam+2./3.*mu)*alpha*(T-T_ref)*kronecker(mydomain)
    mypde.setValue(A=C,X=sigma0,q=msk)
    mypde.getSolverOptions().setVerbosityOn()
    #... solve pde ...
    u=mypde.getSolution()
    #... calculate von-Misses
    g=grad(u)
    sigma=mu*(g+transpose(g))+lam*trace(g)*kronecker(mydomain)-sigma0
    sigma_mises=sqrt(((sigma[0,0]-sigma[1,1])**2+(sigma[1,1]-sigma[2,2])**2+ \
                      (sigma[2,2]-sigma[0,0])**2)/2. \
                       +3*(sigma[0,1]**2 + sigma[1,2]**2 + sigma[2,0]**2))
    #... output ...
    saveVTK("output/deform.vtu",disp=u,stress=sigma_mises)
 
