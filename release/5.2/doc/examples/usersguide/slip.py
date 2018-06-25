from __future__ import division
##############################################################################
#
# Copyright (c) 2009-2018 by The University of Queensland
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

__copyright__="""Copyright (c) 2009-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"


from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.escript.models import FaultSystem
try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
from esys.weipa import saveVTK
from esys.escript.unitsSI import DEG

if not HAVE_FINLEY:
    print("Finley module not available")
else:
    #... set some parameters ...
    lam=1.
    mu=1
    slip_max=1.

    mydomain = Rectangle(l0=1.,l1=1.,n0=16, n1=16)  # n1 need to be multiple of 4!!!
    # .. create the fault system
    fs=FaultSystem(dim=2)
    fs.addFault(V0=[0.5,0.25], strikes=90*DEG, ls=0.5, tag=1)
    # ... create a slip distribution on the fault:
    p, m=fs.getParametrization(mydomain.getX(),tag=1)
    p0,p1= fs.getW0Range(tag=1)
    s=m*(p-p0)*(p1-p)/((p1-p0)/2)**2*slip_max*[0.,1.]
    # ... calculate stress according to slip:
    D=symmetric(grad(s))
    chi, d=fs.getSideAndDistance(D.getFunctionSpace().getX(),tag=1)
    sigma_s=(mu*D+lam*trace(D)*kronecker(mydomain))*chi
    #... open symmetric PDE ...
    mypde=LinearPDE(mydomain)
    mypde.setSymmetryOn()
    #... set coefficients ...
    C=Tensor4(0.,Function(mydomain))
    for i in range(mydomain.getDim()):
      for j in range(mydomain.getDim()):
         C[i,i,j,j]+=lam
         C[j,i,j,i]+=mu
         C[j,i,i,j]+=mu
    # ... fix displacement in normal direction 
    x=mydomain.getX()
    msk=whereZero(x[0])*[1.,0.] + whereZero(x[0]-1.)*[1.,0.] \
       +whereZero(x[1])*[0.,1.] + whereZero(x[1]-1.)*[0.,1.]
    mypde.setValue(A=C,X=-0.5*sigma_s,q=msk)
    #... solve pde ...
    mypde.getSolverOptions().setVerbosityOn()
    v=mypde.getSolution()
    # .. write the displacement to file:
    D=symmetric(grad(v))
    sigma=(mu*D+lam*trace(D)*kronecker(mydomain))+0.5*sigma_s
    saveVTK("slip.vtu",disp=v+0.5*chi*s, stress= sigma)
