
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
import tempfile
      
from esys.escript import *
from esys.finley import Rectangle
from esys.escript.models import DarcyFlow
import sys
import os
try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'


VERBOSE=False # or True
DETAIL_VERBOSE=False

from esys.escript import *
from esys.escript.models import StokesProblemCartesian, PowerLaw
from esys.finley import Rectangle, Brick

from esys.escript.models import Mountains
from math import pi

#====================================================================================================================
class Test_StokesProblemCartesian2D(unittest.TestCase):
   def setUp(self):
       NE=6
       self.TOL=1e-3
       self.domain=Rectangle(NE,NE,order=2,useFullElementOrder=True)
   def tearDown(self):
       del self.domain
   def test_PCG_P_0(self):
       ETA=1.
       P1=0.

       x=self.domain.getX()
       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
       mask=whereZero(x[0])    * [1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.] \
              +whereZero(x[1])    * [1.,0.] \
              +whereZero(x[1]-1)  * [1.,1.]
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*[0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE,max_iter=100,usePCG=True)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       error_p=Lsup(p+P1*x[0]*x[1])
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")

   def test_PCG_P_small(self):
       ETA=1.
       P1=1.

       x=self.domain.getX()
       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
       mask=whereZero(x[0])    * [1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.] \
              +whereZero(x[1])    * [1.,0.] \
              +whereZero(x[1]-1)  * [1.,1.]
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*[0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.2)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE,max_iter=100,usePCG=True)
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       error_p=Lsup(P1*x[0]*x[1]+p)
       saveVTK("d.vtu",p=p, e=P1*x[0]*x[1]+p, p_ref=P1*x[0]*x[1])
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")

   def test_PCG_P_large(self):
       ETA=1.
       P1=1000.

       x=self.domain.getX()
       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
       mask=whereZero(x[0])    * [1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.] \
              +whereZero(x[1])    * [1.,0.] \
              +whereZero(x[1]-1)  * [1.,1.]
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*[0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE,max_iter=100,usePCG=True)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       error_p=Lsup(P1*x[0]*x[1]+p)/P1
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")

   def test_GMRES_P_0(self):
       ETA=1.
       P1=0.

       x=self.domain.getX()
       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
       mask=whereZero(x[0])    * [1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.] \
              +whereZero(x[1])    * [1.,0.] \
              +whereZero(x[1]-1)  * [1.,1.]
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*[0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE,max_iter=50,usePCG=False,iter_restart=18)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       error_p=Lsup(P1*x[0]*x[1]+p)
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")

   def test_GMRES_P_small(self):
       ETA=1.
       P1=1.

       x=self.domain.getX()
       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
       mask=whereZero(x[0])    * [1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.] \
              +whereZero(x[1])    * [1.,0.] \
              +whereZero(x[1]-1)  * [1.,1.]
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*[0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.1)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE,max_iter=20,usePCG=False)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       error_p=Lsup(P1*x[0]*x[1]+p)
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")

   def test_GMRES_P_large(self):
       ETA=1.
       P1=1000.

       x=self.domain.getX()
       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
       mask=whereZero(x[0])    * [1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.] \
              +whereZero(x[1])    * [1.,0.] \
              +whereZero(x[1]-1)  * [1.,1.]
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*[0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE,max_iter=100,usePCG=False)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       error_p=Lsup(P1*x[0]*x[1]+p)/P1
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
#====================================================================================================================
class Test_StokesProblemCartesian3D(unittest.TestCase):
   def setUp(self):
       NE=6
       self.TOL=1e-4
       self.domain=Brick(NE,NE,NE,order=2,useFullElementOrder=True)
   def tearDown(self):
       del self.domain
   def test_PCG_P_0(self):
       ETA=1.
       P1=0.

       x=self.domain.getX()
       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
       x=self.domain.getX()
       mask=whereZero(x[0])    * [1.,1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.,1.] \
              +whereZero(x[1])    * [1.,0.,1.] \
              +whereZero(x[1]-1)  * [1.,1.,1.] \
              +whereZero(x[2])    * [1.,1.,0.] \
              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE ,max_iter=100,usePCG=True)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       error_p=Lsup(P1*x[0]*x[1]*x[2]+p)
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
       self.failUnless(error_v2<10*self.TOL, "2-velocity error too large.")

   def test_PCG_P_small(self):
       ETA=1.
       P1=1.

       x=self.domain.getX()
       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
       mask=whereZero(x[0])    * [1.,1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.,1.] \
              +whereZero(x[1])    * [1.,0.,1.] \
              +whereZero(x[1]-1)  * [1.,1.,1.] \
              +whereZero(x[2])    * [1.,1.,0.] \
              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.1)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE ,max_iter=100,usePCG=True)
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       error_p=Lsup(P1*x[0]*x[1]*x[2]+p)
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
       self.failUnless(error_v2<10*self.TOL, "2-velocity error too large.")
   def test_PCG_P_large(self):
       ETA=1.
       P1=1000.

       x=self.domain.getX()
       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
       mask=whereZero(x[0])    * [1.,1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.,1.] \
              +whereZero(x[1])    * [1.,0.,1.] \
              +whereZero(x[1]-1)  * [1.,1.,1.] \
              +whereZero(x[2])    * [1.,1.,0.] \
              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,-p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE ,max_iter=100,usePCG=True)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       error_p=Lsup(P1*x[0]*x[1]*x[2]+p)/P1 
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
       self.failUnless(error_v2<10*self.TOL, "2-velocity error too large.")

   def test_GMRES_P_0(self):
       ETA=1.
       P1=0.

       x=self.domain.getX()
       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
       x=self.domain.getX()
       mask=whereZero(x[0])    * [1.,1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.,1.] \
              +whereZero(x[1])    * [1.,1.,1.] \
              +whereZero(x[1]-1)  * [1.,1.,1.] \
              +whereZero(x[2])    * [1.,1.,0.] \
              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE,max_iter=100,usePCG=False,iter_restart=20)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       error_p=Lsup(P1*x[0]*x[1]*x[2]+p)
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
       self.failUnless(error_v2<10*self.TOL, "2-velocity error too large.")
   def test_GMRES_P_small(self):
       ETA=1.
       P1=1.

       x=self.domain.getX()
       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
       mask=whereZero(x[0])    * [1.,1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.,1.] \
              +whereZero(x[1])    * [1.,1.,1.] \
              +whereZero(x[1]-1)  * [1.,1.,1.] \
              +whereZero(x[2])    * [1.,1.,0.] \
              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.1)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE and False, verbose=VERBOSE,max_iter=100,usePCG=False)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       error_p=Lsup(P1*x[0]*x[1]*x[2]+p)/P1
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
       self.failUnless(error_v2<10*self.TOL, "2-velocity error too large.")
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
   def test_GMRES_P_large(self):
       ETA=1.
       P1=1000.

       x=self.domain.getX()
       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
       mask=whereZero(x[0])    * [1.,1.,1.] \
              +whereZero(x[0]-1)  * [1.,1.,1.] \
              +whereZero(x[1])    * [1.,0.,1.] \
              +whereZero(x[1]-1)  * [1.,1.,1.] \
              +whereZero(x[2])    * [1.,1.,0.] \
              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
       sp=StokesProblemCartesian(self.domain)
       
       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
       p0=Scalar(P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,show_details=DETAIL_VERBOSE, verbose=VERBOSE ,max_iter=100,usePCG=False)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       error_p=Lsup(P1*x[0]*x[1]*x[2]+p)/P1
       self.failUnless(error_p<10*self.TOL, "pressure error too large.")
       self.failUnless(error_v0<10*self.TOL, "0-velocity error too large.")
       self.failUnless(error_v1<10*self.TOL, "1-velocity error too large.")
       self.failUnless(error_v2<10*self.TOL, "2-velocity error too large.")
#====================================================================================================================
class Test_Darcy(unittest.TestCase):
    # this is a simple test for the darcy flux problem
    #
    # 
    #  p = 1/k * ( 1/2* (fb-f0)/xb* x **2 + f0 * x - ub*x ) +  p0
    # 
    #  with f = (fb-f0)/xb* x + f0 
    #
    #    u = f - k * p,x = ub
    #
    #  we prescribe pressure at x=x0=0 to p0
    # 
    #  if we prescribe the pressure on the bottom x=xb we set
    # 
    #  pb= 1/k * ( 1/2* (fb-f0)* xb + f0 * xb - ub*xb ) +  p0 = 1/k * ((fb+f0)/2 - ub ) * xb  + p0
    # 
    #  which leads to ub = (fb+f0)/2-k*(pb-p0)/xb
    #
    def rescaleDomain(self):
        x=self.dom.getX().copy()
        for i in xrange(self.dom.getDim()):
             x_inf=inf(x[i])
             x_sup=sup(x[i])
             if i == self.dom.getDim()-1:
                x[i]=-self.WIDTH*(x[i]-x_sup)/(x_inf-x_sup)
             else:
                x[i]=self.WIDTH*(x[i]-x_inf)/(x_sup-x_inf)
        self.dom.setX(x)
    def getScalarMask(self,include_bottom=True):
        x=self.dom.getX().copy()
        x_inf=inf(x[self.dom.getDim()-1])
        x_sup=sup(x[self.dom.getDim()-1])
        out=whereZero(x[self.dom.getDim()-1]-x_sup)
        if include_bottom: out+=whereZero(x[self.dom.getDim()-1]-x_inf)
        return wherePositive(out)
    def getVectorMask(self,include_bottom=True):
        x=self.dom.getX().copy()
        out=Vector(0.,Solution(self.dom))
        for i in xrange(self.dom.getDim()):
             x_inf=inf(x[i])
             x_sup=sup(x[i])
             if i != self.dom.getDim()-1: out[i]+=whereZero(x[i]-x_sup)
             if i != self.dom.getDim()-1 or include_bottom: out[i]+=whereZero(x[i]-x_inf)
        return wherePositive(out)

    def setSolutionFixedBottom(self, p0, pb, f0, fb, k):
         d=self.dom.getDim()
         x=self.dom.getX()[d-1]
         xb=inf(x)
         u=Vector(0.,Solution(self.dom))+kronecker(d)[d-1]*((f0+fb)/2.-k*(pb-p0)/xb)
         p=1./k*((fb-f0)/(xb*2.)* x**2 - (fb-f0)/2.*x)+(pb-p0)/xb*x +  p0
         f= ((fb-f0)/xb* x + f0)*kronecker(Function(self.dom))[d-1]
         return u,p,f
        
    def testConstF_FixedBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u_ref,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")
    def testConstF_FixedBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p,max_iter=100, verbose=VERBOSE )
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")

    def testConstF_FixedBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FixedBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FixedBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FixedBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=True)
        mv=self.getVectorMask(include_bottom=False)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                    location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testConstF_FreeBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=10.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FreeBottom_smallK(self):
        k=1.e-10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*25.*Lsup(u_ref), "flux error too big.")  # 25 because of disc. error.
        self.failUnless(Lsup(p-p_ref)<self.TOL*25.*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FreeBottom_mediumK(self):
        k=1.
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

    def testVarioF_FreeBottom_largeK(self):
        k=1.e10
        mp=self.getScalarMask(include_bottom=False)
        mv=self.getVectorMask(include_bottom=True)
        u_ref,p_ref,f=self.setSolutionFixedBottom(p0=2500,pb=4000.,f0=10.,fb=30.,k=k)
        p=p_ref*mp
        u=u_ref*mv
        df=DarcyFlow(self.dom)
        df.setValue(g=f,
                      location_of_fixed_pressure=mp,
                      location_of_fixed_flux=mv,
                      permeability=Scalar(k,Function(self.dom)))
        df.setTolerance(rtol=self.TOL)
        df.setSubProblemTolerance()
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE)
        self.failUnless(Lsup(v-u_ref)<self.TOL*10.*Lsup(u_ref), "flux error too big.")
        self.failUnless(Lsup(p-p_ref)<self.TOL*10.*Lsup(p_ref), "pressure error too big.")

class Test_Darcy2D(Test_Darcy):
    TOL=1e-4
    WIDTH=1.
    def setUp(self):
        NE=40  # wrning smaller NE may case a failure for VarioF tests due to discretization errors.
        self.dom = Rectangle(NE,NE)
        self.rescaleDomain()
    def tearDown(self):
        del self.dom
class Test_Darcy3D(Test_Darcy):
    TOL=1e-4
    WIDTH=1.
    def setUp(self):
        NE=25  # wrning smaller NE may case a failure for VarioF tests due to discretization errors.
        self.dom = Brick(NE,NE,NE)
        self.rescaleDomain()
    def tearDown(self):
        del self.dom


class Test_Mountains3D(unittest.TestCase):
   def setUp(self):
       NE=16
       self.TOL=1e-4
       self.domain=Brick(NE,NE,NE,order=1,useFullElementOrder=True)
   def tearDown(self):
       del self.domain
   def test_periodic(self):
       EPS=0.01

       x=self.domain.getX()
       v = Vector(0.0, Solution(self.domain))
       a0=1
       a1=1
       n0=2
       n1=2
       n2=0.5
       a2=-(a0*n0+a1*n1)/n2
       v[0]=a0*sin(pi*n0*x[0])* cos(pi*n1*x[1])* cos(pi*n2*x[2])
       v[1]=a1*cos(pi*n0*x[0])* sin(pi*n1*x[1])* cos(pi*n2*x[2])
       v[2]=a2*cos(pi*n0*x[0])* cos(pi*n1*x[1])* sin(pi*n2*x[2])

       H_t=Scalar(0.0, Solution(self.domain))
       mts=Mountains(self.domain,v,eps=EPS,z=1)
       u,Z=mts.update(u=v,H_t=H_t)
       
       error_int=integrate(Z)
       self.failUnless(error_int<self.TOL, "Boundary intergral is too large.")

class Test_Mountains2D(unittest.TestCase):
   def setUp(self):
       NE=16
       self.TOL=1e-4
       self.domain=Rectangle(NE,NE,order=1,useFullElementOrder=True)
   def tearDown(self):
       del self.domain
   def test_periodic(self):
       EPS=0.01

       x=self.domain.getX()
       v = Vector(0.0, Solution(self.domain))
       a0=1
       n0=1
       n1=0.5
       a1=-(a0*n0)/n1
       v[0]=a0*sin(pi*n0*x[0])* cos(pi*n1*x[1])
       v[1]=a1*cos(pi*n0*x[0])* sin(pi*n1*x[1])
       
       H_t=Scalar(0.0, Solution(self.domain))
       mts=Mountains(self.domain,v,eps=EPS,z=1)
       u,Z=mts.update(u=v,H_t=H_t)
       
       error_int=integrate(Z)
       self.failUnless(error_int<self.TOL, "Boundary intergral is too large.")
       


class Test_Rheologies(unittest.TestCase):

     """
     this is the program used to generate the powerlaw tests:

     TAU_Y=100.
     N=10
     M=5

     def getE(tau):
         if tau<=TAU_Y:
           return 1./(0.5+20*sqrt(tau))
         else:
           raise ValueError,"out of range."
     tau=[ i* (TAU_Y/N) for i in range(N+1)] + [TAU_Y for i in range(M)]
     e=[ tau[i]/getE(tau[i]) for i in range(N+1)]
     e+= [ (i+1+j)* (max(e)/N) for j in range(M) ]

     print tau
     print e
     """
     TOL=1.e-8
     def test_PowerLaw_Init(self):
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)

         self.failUnlessEqual(pl.getNumMaterials(),3,"num materials is wrong")
         self.failUnlessEqual(pl.validMaterialId(0),True,"material id 0 not found")
         self.failUnlessEqual(pl.validMaterialId(1),True,"material id 1 not found")
         self.failUnlessEqual(pl.validMaterialId(2),True,"material id 2 not found")
         self.failUnlessEqual(pl.validMaterialId(3),False,"material id 3 not found")

         self.failUnlessEqual(pl.getFriction(),None,"initial friction wrong.")
         self.failUnlessEqual(pl.getTauY(),None,"initial tau y wrong.")
         pl.setDruckerPragerLaw(tau_Y=10,friction=3)
         self.failUnlessEqual(pl.getFriction(),3,"friction wrong.")
         self.failUnlessEqual(pl.getTauY(),10,"tau y wrong.")

         self.failUnlessEqual(pl.getElasticShearModulus(),None,"initial shear modulus is wrong")
         pl.setElasticShearModulus(1000)
         self.failUnlessEqual(pl.getElasticShearModulus(),1000.,"shear modulus is wrong")

         e=pl.getEtaN()
         self.failUnlessEqual(len(e),3,"initial length of etaN is wrong.")
         self.failUnlessEqual(e,[None, None, None],"initial etaN are wrong.")
         self.failUnlessEqual(pl.getEtaN(0),None,"initial material id 0 is wrong")
         self.failUnlessEqual(pl.getEtaN(1),None,"initial material id 1 is wrong")
         self.failUnlessEqual(pl.getEtaN(2),None,"initial material id 2 is wrong")
         self.failUnlessRaises(ValueError, pl.getEtaN, 3)

         self.failUnlessRaises(ValueError,pl.setPowerLaws,eta_N=[10,20,30,40],tau_t=[2], power=[5,6,7])
         self.failUnlessRaises(ValueError,pl.setPowerLaws,eta_N=[10,20,30],tau_t=[2,4,5,5], power=[5,6,7])
         self.failUnlessRaises(ValueError,pl.setPowerLaws,eta_N=[10,20,30],tau_t=[2,1,2], power=[5,6,7,7])

         pl.setPowerLaws(eta_N=[10,20,30],tau_t=[100,200,300], power=[1,2,3])
         self.failUnlessEqual(pl.getPower(),[1,2,3],"powers are wrong.")
         self.failUnlessEqual(pl.getTauT(),[100,200,300],"tau t are wrong.")
         self.failUnlessEqual(pl.getEtaN(),[10,20,30],"etaN are wrong.")

         pl.setPowerLaw(id=0,eta_N=40,tau_t=400, power=4)
         self.failUnlessEqual(pl.getPower(),[4,2,3],"powers are wrong.")
         self.failUnlessEqual(pl.getTauT(),[400,200,300],"tau t are wrong.")
         self.failUnlessEqual(pl.getEtaN(),[40,20,30],"etaN are wrong.")

         self.failUnlessRaises(ValueError,pl.getPower,-1)
         self.failUnlessEqual(pl.getPower(0),4,"power 0 is wrong.")
         self.failUnlessEqual(pl.getPower(1),2,"power 1 is wrong.")
         self.failUnlessEqual(pl.getPower(2),3,"power 2 is wrong.")
         self.failUnlessRaises(ValueError,pl.getPower,3)

         self.failUnlessRaises(ValueError,pl.getTauT,-1)
         self.failUnlessEqual(pl.getTauT(0),400,"tau t 0 is wrong.")
         self.failUnlessEqual(pl.getTauT(1),200,"tau t 1 is wrong.")
         self.failUnlessEqual(pl.getTauT(2),300,"tau t 2 is wrong.")
         self.failUnlessRaises(ValueError,pl.getTauT,3)

         self.failUnlessRaises(ValueError,pl.getEtaN,-1)
         self.failUnlessEqual(pl.getEtaN(0),40,"eta n 0 is wrong.")
         self.failUnlessEqual(pl.getEtaN(1),20,"eta n 1 is wrong.")
         self.failUnlessEqual(pl.getEtaN(2),30,"eta n 2 is wrong.")
         self.failUnlessRaises(ValueError,pl.getEtaN,3)

     def checkResult(self,id,gamma_dot_, eta, tau_ref):
         self.failUnless(eta>=0,"eta needs to be positive (test %s)"%id)
         error=abs(gamma_dot_*eta-tau_ref)
         self.failUnless(error<=self.TOL*tau_ref,"eta is wrong: error = gamma_dot_*eta-tau_ref = %s * %s - %s = %s (test %s)"%(gamma_dot_,eta,tau_ref,error,id))
        
     def test_PowerLaw_Linear(self):
         taus= [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0]
         pl=PowerLaw(numMaterials=1,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaw(eta_N=2.)
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])
        
     def test_PowerLaw_QuadLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 637.45553203367592, 1798.8543819998317, 3301.3353450309969, 5079.6442562694074, 7096.067811865476, 9325.1600308977995, 11748.240371477059, 14350.835055998656, 17121.29936490925, 20050.0, 22055.0, 24060.0, 26065.0, 28070.0, 30075.0]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01],tau_t=[1, 25.], power=[1,2])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.632455532033676, 11.788854381999831, 18.286335345030995, 25.059644256269408, 32.071067811865476, 39.295160030897804, 46.713240371477056, 54.310835055998652, 62.076299364909254, 70.0, 77.0, 84.0, 91.0, 98.0, 105.0]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.],tau_t=[1, 25.], power=[1,2])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 51.415888336127786, 157.36125994561547, 304.6468153816889, 487.84283811405851, 703.60440414872664, 949.57131887226353, 1223.9494765692673, 1525.3084267560891, 1852.4689652218574, 2204.4346900318833, 2424.8781590350718, 2645.3216280382603, 2865.7650970414484, 3086.2085660446369, 3306.6520350478249]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,1./16.],tau_t=[1, 64.], power=[1,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.4641588833612778, 11.473612599456157, 17.89646815381689, 24.678428381140588, 31.786044041487269, 39.195713188722635, 46.889494765692675, 54.853084267560895, 63.074689652218574, 71.544346900318828, 78.698781590350706, 85.853216280382583, 93.007650970414474, 100.16208566044635, 107.316520350478]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,25./4.],tau_t=[1, 64.], power=[1,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadLarge_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 683.87142036980367, 1946.2156419454475, 3590.982160412686, 5547.4870943834667, 7774.6722160142008, 10244.731349770063, 12937.189848046326, 15836.143482754744, 18928.768330131104, 22204.434690031885, 24424.878159035074, 26645.321628038262, 28865.765097041451, 31086.208566044639, 33306.652035047824]
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01,1./16.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadLarge_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 637.9196909170372, 1800.3279945992881, 3304.2318131848137, 5084.3226846505486, 7102.853855906963, 9334.3557440865225, 11760.129866242751, 14365.688140266215, 17139.374054561471, 20071.544346900318, 22078.698781590349, 24085.853216280382, 26093.007650970416, 28100.162085660446, 30107.316520350476]
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01,25./4.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 52.04834386816146, 159.15011432761528, 307.93315072671987, 492.9024823703279, 710.67547196059206, 958.86647890316135, 1235.6627169407443, 1539.6192618120876, 1869.5452645867665, 2224.4346900318833, 2446.8781590350718, 2669.3216280382603, 2891.7650970414484, 3114.2085660446369, 3336.6520350478249]

         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.,1./16.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 6.0966144153949529, 13.262466981455987, 21.182803498847885, 29.738072637409996, 38.857111853352741, 48.49087321962044, 58.602735137169738, 69.16391932355954, 80.150989017127827, 91.544346900318828, 100.69878159035071, 109.85321628038258, 119.00765097041447, 128.16208566044637, 137.31652035047824]
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.,25./4.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_withShear(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0]
         pl=PowerLaw(numMaterials=1,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaw(eta_N=2.)
         pl.setElasticShearModulus(3.)
         dt=1./3.
         pl.setEtaTolerance(self.TOL)
         self.failUnlessRaises(ValueError, pl.getEtaEff,gamma_dot_s[0])
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i],dt=dt),taus[i])

if __name__ == '__main__':
   suite = unittest.TestSuite()
   
   suite.addTest(unittest.makeSuite(Test_StokesProblemCartesian2D))
   suite.addTest(unittest.makeSuite(Test_Darcy3D))
   suite.addTest(unittest.makeSuite(Test_Darcy2D))
   suite.addTest(unittest.makeSuite(Test_StokesProblemCartesian3D))
   suite.addTest(unittest.makeSuite(Test_Mountains3D))
   suite.addTest(unittest.makeSuite(Test_Mountains2D))
   suite.addTest(unittest.makeSuite(Test_Rheologies))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

