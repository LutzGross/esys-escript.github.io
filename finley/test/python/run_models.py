
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
__url__="http://www.uq.edu.au/esscc/escript-finley"

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


VERBOSE=False  # or True
DETAIL_VERBOSE=False

from esys.escript import *
from esys.escript.models import StokesProblemCartesian
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u_ref,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p,max_iter=100, verbose=VERBOSE ,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
        df.setTolerance(rtol=self.TOL,p_ref=p_ref,v_ref=u_ref)
        v,p=df.solve(u,p, max_iter=100, verbose=VERBOSE,sub_rtol=self.TOL/200)
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
       
if __name__ == '__main__':
   suite = unittest.TestSuite()
   
   suite.addTest(unittest.makeSuite(Test_StokesProblemCartesian2D))
   suite.addTest(unittest.makeSuite(Test_Darcy3D))
   suite.addTest(unittest.makeSuite(Test_Darcy2D))
   suite.addTest(unittest.makeSuite(Test_StokesProblemCartesian3D))
   suite.addTest(unittest.makeSuite(Test_Mountains3D))
   suite.addTest(unittest.makeSuite(Test_Mountains2D))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

