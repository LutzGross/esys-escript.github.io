
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
import tempfile
      


VERBOSE=False # or True

from esys.escript import *
from esys.escript.models import StokesProblemCartesian, PowerLaw, IncompressibleIsotropicFlowCartesian, FaultSystem, DarcyFlow
from esys.escript.models import Mountains
from esys.finley import Rectangle, Brick

from math import pi
import numpy
import sys
import os
#====================================================================================================================
try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0,verbose=VERBOSE,max_iter=100,usePCG=True)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.2)
       u,p=sp.solve(u0,p0, verbose=VERBOSE,max_iter=100,usePCG=True)
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       error_p=Lsup(P1*x[0]*x[1]+p)
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0, verbose=VERBOSE,max_iter=100,usePCG=True)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0, verbose=VERBOSE,max_iter=50,usePCG=False,iter_restart=18)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.1)
       u,p=sp.solve(u0,p0, verbose=VERBOSE,max_iter=20,usePCG=False)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0, verbose=VERBOSE,max_iter=100,usePCG=False)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0, verbose=VERBOSE ,max_iter=100,usePCG=True)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.1)
       u,p=sp.solve(u0,p0, verbose=VERBOSE ,max_iter=100,usePCG=True)
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,-p0, verbose=VERBOSE ,max_iter=100,usePCG=True)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0, verbose=VERBOSE,max_iter=100,usePCG=False,iter_restart=20)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL*0.1)
       u,p=sp.solve(u0,p0, verbose=VERBOSE,max_iter=100,usePCG=False)
       
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
       p0=Scalar(-P1,ReducedSolution(self.domain))
       sp.setTolerance(self.TOL)
       u,p=sp.solve(u0,p0, verbose=VERBOSE ,max_iter=100,usePCG=False)
       
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

       mts=Mountains(self.domain,eps=EPS)
       mts.setVelocity(v)
       Z=mts.update()
       
       error_int=abs(integrate(Z*whereZero(FunctionOnBoundary(self.domain).getX()[2]-1.)))
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
       mts=Mountains(self.domain,eps=EPS)
       mts.setVelocity(v)
       Z=mts.update()
       
       error_int=abs(integrate(Z*whereZero(FunctionOnBoundary(self.domain).getX()[1]-1.)))
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
         gamma_dot_s=[0.0, 405.0, 1610.0, 3615.0, 6420.0, 10025.0, 14430.0, 19635.0, 25640.0, 32445.0, 40050.0, 44055.0, 48060.0, 52065.0, 56070.0, 60075.0]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01],tau_t=[1, 25.], power=[1,2])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.4, 11.6, 18.6, 26.4, 35.0, 44.4, 54.6, 65.6, 77.4, 90.0, 99.0, 108.0, 117.0, 126.0, 135.0]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.],tau_t=[1, 25.], power=[1,2])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 8.90625, 41.25, 120.46875, 270.0, 513.28125, 873.75, 1374.84375, 2040.0, 2892.65625, 3956.25, 4351.875, 4747.5, 5143.125, 5538.75, 5934.375]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,1./16.],tau_t=[1, 64.], power=[1,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.0390625, 10.3125, 16.0546875, 22.5, 29.8828125, 38.4375, 48.3984375, 60.0, 73.4765625, 89.0625, 97.96875, 106.875, 115.78125, 124.6875, 133.59375]
         pl=PowerLaw(numMaterials=2,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,25./4.],tau_t=[1, 64.], power=[1,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadLarge_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 408.90625, 1641.25, 3720.46875, 6670.0, 10513.28125, 15273.75, 20974.84375, 27640.000000000004, 35292.65625, 43956.25, 48351.875, 52747.5, 57143.125, 61538.75, 65934.375]

         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01,1./16.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadLarge_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 405.0390625, 1610.3125, 3616.0546875, 6422.5, 10029.8828125, 14438.4375, 19648.3984375, 25660.0, 32473.4765625, 40089.0625, 44097.96875, 48106.875, 52115.78125, 56124.6875, 60133.59375]

         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,0.01,25./4.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall_CubeLarge(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 9.30625, 42.85, 124.06875, 276.4, 523.28125, 888.15, 1394.44375, 2065.6, 2925.05625, 3996.25, 4395.875, 4795.5, 5195.125, 5594.75, 5994.375]
         pl=PowerLaw(numMaterials=3,verbose=VERBOSE)
         pl.setDruckerPragerLaw(tau_Y=100.)
         pl.setPowerLaws(eta_N=[2.,10.,1./16.],tau_t=[1, 25.,64.], power=[1,2,3])
         pl.setEtaTolerance(self.TOL)
         for i in xrange(len(taus)): self.checkResult(i,gamma_dot_s[i], pl.getEtaEff(gamma_dot_s[i]),taus[i])

     def test_PowerLaw_QuadSmall_CubeSmall(self):
         taus=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         gamma_dot_s=[0.0, 5.4390625, 11.9125, 19.6546875, 28.9, 39.8828125, 52.8375, 67.9984375, 85.6, 105.8765625, 129.0625, 141.96875, 154.875, 167.78125, 180.6875, 193.59375]
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

class Test_IncompressibleIsotropicFlowCartesian(unittest.TestCase):
   TOL=1.e-6
   VERBOSE=False
   A=1.
   P_max=100
   NE=2*getMPISizeWorld()
   tau_Y=10.
   N_dt=10

   # material parameter:
   tau_1=5.
   tau_2=5.
   eta_0=100.
   eta_1=50.
   eta_2=400.
   N_1=2.
   N_2=3.
   def getReference(self, t):

      B=self.tau_Y/sqrt((self.dom.getDim()-1)*self.dom.getDim()*0.5)
      x=self.dom.getX()

      s_00=min(self.A*t,B)
      tau=sqrt((self.dom.getDim()-1)*self.dom.getDim()*0.5)*abs(s_00)
      inv_eta= 1./self.eta_0 + 1./self.eta_1*(tau/self.tau_1)**(self.N_1-1.) + 1./self.eta_2*(tau/self.tau_2)**(self.N_2-1.)

      alpha=0.5*inv_eta*s_00
      if s_00 <= B and self.mu !=None: alpha+=1./(2*self.mu)*self.A
      u_ref=x*alpha
      u_ref[self.dom.getDim()-1]=(1.-x[self.dom.getDim()-1])*alpha*(self.dom.getDim()-1)
      sigma_ref=kronecker(self.dom)*s_00
      sigma_ref[self.dom.getDim()-1,self.dom.getDim()-1]=-s_00*(self.dom.getDim()-1)

      p_ref=self.P_max
      for d in range(self.dom.getDim()): p_ref=p_ref*x[d]
      p_ref-=integrate(p_ref)/vol(self.dom)
      return u_ref, sigma_ref, p_ref

   def runIt(self, free=None):
      x=self.dom.getX()
      B=self.tau_Y/sqrt((self.dom.getDim()-1)*self.dom.getDim()*0.5)
      dt=B/int(self.N_dt/2)
      if self.VERBOSE: print "dt =",dt
      if self.latestart:
          t=dt
      else:
          t=0
      v,s,p=self.getReference(t)

      mod=IncompressibleIsotropicFlowCartesian(self.dom, stress=s, v=v, p=p, t=t, numMaterials=3, verbose=self.VERBOSE)
      mod.setDruckerPragerLaw(tau_Y=self.tau_Y,friction=None)
      mod.setElasticShearModulus(self.mu)
      mod.setPowerLaws([self.eta_0, self.eta_1, self.eta_2], [ 1., self.tau_1, self.tau_2],  [1.,self.N_1,self.N_2])
      mod.setTolerance(self.TOL)
      mod.setFlowTolerance(self.TOL)

      BF=Vector(self.P_max,Function(self.dom))
      for d in range(self.dom.getDim()):
          for d2 in range(self.dom.getDim()):
              if d!=d2: BF[d]*=x[d2]
      v_mask=Vector(0,Solution(self.dom))
      if free==None:
         for d in range(self.dom.getDim()):
            v_mask[d]=whereZero(x[d])+whereZero(x[d]-1.)
      else:
         for d in range(self.dom.getDim()):
            if d == self.dom.getDim()-1:
               v_mask[d]=whereZero(x[d]-1.)
            else:
               v_mask[d]=whereZero(x[d])
      mod.setExternals(F=BF,fixed_v_mask=v_mask)
       
      n=self.dom.getNormal()
      N_t=0
      errors=[]
      while N_t < self.N_dt:
         t_ref=t+dt
         v_ref, s_ref,p_ref=self.getReference(t_ref)
         mod.setExternals(f=matrixmult(s_ref,n)-p_ref*n, v_boundary=v_ref)
         mod.update(dt, iter_max=100, inner_iter_max=20, verbose=self.VERBOSE, usePCG=True)
         self.check(N_t,mod,t_ref,v_ref, s_ref,p_ref)
         t+=dt
         N_t+=1

   def check(self,N_t,mod,t_ref,v_ref, s_ref,p_ref):
         p=mod.getPressure()
         p-=integrate(p)/vol(self.dom)
         error_p=Lsup(mod.getPressure()-p_ref)/Lsup(p_ref)
         error_s=Lsup(mod.getDeviatoricStress()-s_ref)/Lsup(s_ref)
         error_v=Lsup(mod.getVelocity()-v_ref)/Lsup(v_ref)
         error_t=abs(mod.getTime()-t_ref)/abs(t_ref)
         if self.VERBOSE: print "time step ",N_t,"time = ",mod.getTime(),"errors s,p,v = ",error_s, error_p, error_v
         self.failUnless( error_p <= 10*self.TOL, "time step %s: pressure error %s too high."%(N_t,error_p) )
         self.failUnless( error_s <= 10*self.TOL, "time step %s: stress error %s too high."%(N_t,error_s) )
         self.failUnless( error_v <= 10*self.TOL, "time step %s: velocity error %s too high."%(N_t,error_v) )
         self.failUnless( error_t <= 10*self.TOL, "time step %s: time marker error %s too high."%(N_t,error_t) )
   def tearDown(self):
        del self.dom

   def test_D2_Fixed_MuNone_LateStart(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=True
       self.runIt()
   def test_D2_Fixed_Mu_LateStart(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=True
       self.runIt()
   def test_D2_Fixed_MuNone(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=False
       self.runIt()
   def test_D2_Fixed_Mu(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=False
       self.runIt()
   def test_D2_Free_MuNone_LateStart(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=True
       self.runIt(free=0)
   def test_D2_Free_Mu_LateStart(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=True
       self.runIt(free=0)
   def test_D2_Free_MuNone(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=False
       self.runIt(free=0)
   def test_D2_Free_Mu(self):
       self.dom = Rectangle(self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=False
       self.runIt(free=0)

   def test_D3_Fixed_MuNone_LateStart(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=True
       self.runIt()
   def test_D3_Fixed_Mu_LateStart(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=True
       self.runIt()
   def test_D3_Fixed_MuNone(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=False
       self.runIt()
   def test_D3_Fixed_Mu(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=False
       self.runIt()
   def test_D3_Free_MuNone_LateStart(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=True
       self.runIt(free=0)
   def test_D3_Free_Mu_LateStart(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=True
       self.runIt(free=0)
   def test_D3_Free_MuNone(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=None
       self.latestart=False
       self.runIt(free=0)
   def test_D3_Free_Mu(self):
       self.dom = Brick(self.NE,self.NE,self.NE,order=2)
       self.mu=555.
       self.latestart=False
       self.runIt(free=0)


class Test_FaultSystem(unittest.TestCase):
   EPS=1.e-8
   NE=10
   def test_Fault2D_MaxValue(self):
      dom=Rectangle(2*self.NE,2*self.NE)
      x=dom.getX()
      f=FaultSystem(dim=2)
      f.addFault([[0.5,0.],[0.,0.5]], tag=1)
      f.addFault([[1.,0.5],[0.5,0.5],[0.5,1.]], tag=2)

      m, t, l=f.getMaxValue(x[0]*(1.-x[0])*(1-x[1]))
      self.failUnless(  m == 0.25, "wrong max value")
      self.failUnless(  t == 1, "wrong max tag")
      self.failUnless(  l == 0., "wrong max location")
      m, t, l=f.getMaxValue(x[1]*(1.-x[1])*(1-x[0])*x[0])
      self.failUnless(  m == 0.0625, "wrong max value")
      self.failUnless(  t == 2, "wrong max tag")
      self.failUnless(  l == 0.5, "wrong max location")
      m, t, l=f.getMaxValue(x[0]*(1.-x[0])*x[1])
      self.failUnless(  m == 0.25, "wrong max value")
      self.failUnless(  t == 2, "wrong max tag")
      self.failUnless(  l == 1.0, "wrong max location")
      m, t, l= f.getMaxValue(x[1]*(1.-x[1])*x[0])
      self.failUnless(  m == 0.25, "wrong max value")
      self.failUnless(  t == 2, "wrong max tag")
      self.failUnless(  l == 0., "wrong max location")
      m, t, l= f.getMaxValue(x[1]*(1.-x[1])*(1.-x[0]))
      self.failUnless(  m == 0.25, "wrong max value")
      self.failUnless(  t == 1, "wrong max tag")
      self.failUnless(  abs(l-0.70710678118654) <= self.EPS,  "wrong max location")

      
   def xtest_Fault2D_TwoFaults(self):
      f=FaultSystem(dim=2)
      top1=[ [1.,0.], [1.,1.], [0.,1.] ]
      self.failUnlessRaises(ValueError,f.addFault,top=top1,tag=1,bottom=top1)
      f.addFault(top=top1,tag=1)
      self.failUnless(f.getDim() == 2, "wrong dimension")
      self.failUnless( [ 1 ] == f.getTags(), "tags wrong")
      self.failUnless(  2. == f.getLength(1), "length wrong")
      self.failUnless(  0. == f.getDepth(1), "depth wrong")
      self.failUnless( (0., 2.) ==  f.getW0Range(1)," wrong W0 range")
      self.failUnless( (0., 0.) ==  f.getW1Range(1)," wrong W1 range")
      self.failUnless( [0., 1., 2.] ==  f.getW0Offsets(1)," wrong W0 offsets")
      segs=f.getFaultSegments(1)[0]
      self.failUnless( len(segs) == 3, "wrong number of segments")
      self.failUnless( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.failUnless( segs[0].size == 2, "seg 0 has wrong size.")
      self.failUnless( numpy.linalg.norm(segs[0]-[1.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.failUnless( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.failUnless( segs[1].size == 2, "seg 1 has wrong size.")
      self.failUnless( numpy.linalg.norm(segs[1]-[1.,1.]) < self.EPS, "wrong vertex. 1 ")
      self.failUnless( isinstance(segs[2], numpy.ndarray), "wrong class of vertex 2")
      self.failUnless( segs[2].size == 2, "seg 2 has wrong size.")
      self.failUnless( numpy.linalg.norm(segs[2]-[0.,1.]) < self.EPS, "wrong vertex. 2 ")
      c=f.getCenterOnSurface()
      self.failUnless( isinstance(c, numpy.ndarray), "center has wrong class")
      self.failUnless( c.size == 2, "center size is wrong")
      self.failUnless( numpy.linalg.norm(c-[2./3.,2./3.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.failUnless( abs(o-45.) < self.EPS, "wrong orientation.")

      top2=[ [10.,0.], [0.,10.] ]
      f.addFault(top2, tag=2, w0_offsets=[0,20], w1_max=20)
      self.failUnless( [ 1, 2 ] == f.getTags(), "tags wrong")
      self.failUnless(  abs(f.getLength(2)-14.1421356237) < self.EPS * 14.1421356237, "wrong length")
      self.failUnless(  0. == f.getDepth(2), "depth wrong")
      self.failUnless( (0., 20.) ==  f.getW0Range(2)," wrong W0 range")
      self.failUnless( (0., 0.) ==  f.getW1Range(2)," wrong W1 range")
      self.failUnless( [0., 20.] ==  f.getW0Offsets(2)," wrong W0 offsets")
      segs=f.getFaultSegments(2)[0]
      self.failUnless( len(segs) == 2, "wrong number of segments")
      self.failUnless( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0")
      self.failUnless( numpy.linalg.norm(segs[0]-[10.,0.]) < self.EPS, "wrong vertex. 0 ")
      self.failUnless( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1")
      self.failUnless( numpy.linalg.norm(segs[1]-[0.,10.]) < self.EPS, "wrong vertex. 1 ")
      c=f.getCenterOnSurface()
      self.failUnless( isinstance(c, numpy.ndarray), "center has wrong class")
      self.failUnless( c.size == 2, "center size is wrong")
      self.failUnless( numpy.linalg.norm(c-[12./5.,12./5.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.failUnless( abs(o-45.) < self.EPS, "wrong orientation.")

      f.transform(rot=-pi/2., shift=[-1.,-1.])
      self.failUnless( [ 1, 2 ] == f.getTags(), "tags after transformation wrong")
      self.failUnless(  2. == f.getLength(1), "length after transformation wrong")
      self.failUnless(  0. == f.getDepth(1), "depth after transformation wrong")
      self.failUnless( (0., 2.) ==  f.getW0Range(1)," wrong W0 after transformation range")
      self.failUnless( (0., 0.) ==  f.getW1Range(1)," wrong W1 rangeafter transformation ")
      self.failUnless( [0., 1., 2.] ==  f.getW0Offsets(1)," wrong W0 offsetsafter transformation ")
      segs=f.getFaultSegments(1)[0]
      self.failUnless( len(segs) == 3, "wrong number of segmentsafter transformation ")
      self.failUnless( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0 after transformation")
      self.failUnless( segs[0].size == 2, "seg 0 has wrong size after transformation.")
      self.failUnless( numpy.linalg.norm(segs[0]-[-1.,0.]) < self.EPS, "wrong vertex. 0  after transformation")
      self.failUnless( isinstance(segs[1], numpy.ndarray), "wrong class of vertex  after transformation1")
      self.failUnless( segs[1].size == 2, "seg 1 has wrong size after transformation.")
      self.failUnless( numpy.linalg.norm(segs[1]-[0.,0.]) < self.EPS, "wrong vertex.  after transformation1 ")
      self.failUnless( isinstance(segs[2], numpy.ndarray), "wrong class of vertex  after transformation2")
      self.failUnless( segs[2].size == 2, "seg 2 has wrong size after transformation.")
      self.failUnless( numpy.linalg.norm(segs[2]-[0., 1.]) < self.EPS, "wrong vertex after transformation. 2 ")
      self.failUnless(  abs(f.getLength(2)-14.1421356237) < self.EPS * 14.1421356237, "wrong length after transformation")
      self.failUnless(  0. == f.getDepth(2), "depth wrong after transformation")
      self.failUnless( (0., 20.) ==  f.getW0Range(2)," wrong W0 range after transformation")
      self.failUnless( (0., 0.) ==  f.getW1Range(2)," wrong W1 range after transformation")
      self.failUnless( [0., 20.] ==  f.getW0Offsets(2)," wrong W0 offsets after transformation")
      segs=f.getFaultSegments(2)[0]
      self.failUnless( len(segs) == 2, "wrong number of segments after transformation")
      self.failUnless( isinstance(segs[0], numpy.ndarray), "wrong class of vertex 0 after transformation")
      self.failUnless( numpy.linalg.norm(segs[0]-[-1.,-9]) < self.EPS, "wrong vertex. 0  after transformation")
      self.failUnless( isinstance(segs[1], numpy.ndarray), "wrong class of vertex 1 after transformation")
      self.failUnless( numpy.linalg.norm(segs[1]-[9.,1.]) < self.EPS, "wrong vertex. 1  after transformation")

      c=f.getCenterOnSurface()
      self.failUnless( isinstance(c, numpy.ndarray), "center has wrong class")
      self.failUnless( c.size == 2, "center size is wrong")
      self.failUnless( numpy.linalg.norm(c-[7./5.,-7./5.]) < self.EPS, "center has wrong coordinates.")
      o=f.getOrientationOnSurface()/pi*180.
      self.failUnless( abs(o+45.) < self.EPS, "wrong orientation.")

      p=f.getParametrization([-1.,0.],1)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0])<self.EPS, "wrong value.")
      p=f.getParametrization([-0.5,0.],1)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0]-0.5)<self.EPS* 0.5, "wrong value.")
      p=f.getParametrization([0.,0.],1)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0]-1.)<self.EPS, "wrong value.")
      p=f.getParametrization([0.0000001,0.0000001],1, tol=1.e-8)
      self.failUnless(p[1]==0., "wrong value.")
      p=f.getParametrization([0.0000001,0.0000001],1, tol=1.e-6)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0]-1.0000001)<self.EPS, "wrong value.")
      p=f.getParametrization([0.,0.5],1)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0]-1.5)<self.EPS, "wrong value.")
      p=f.getParametrization([0,1.],1)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0]-2.)<self.EPS, "wrong value.")
      p=f.getParametrization([1.,1.],1)
      self.failUnless(p[1]==0., "wrong value.")
      p=f.getParametrization([0,1.11],1)
      self.failUnless(p[1]==0., "wrong value.")
      p=f.getParametrization([-1,-9.],2)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0])<self.EPS, "wrong value.")
      p=f.getParametrization([9,1],2)
      self.failUnless(p[1]==1., "wrong value.")
      self.failUnless(abs(p[0]-20.)<self.EPS, "wrong value.")




if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_FaultSystem))
   suite.addTest(unittest.makeSuite(Test_StokesProblemCartesian2D))
   suite.addTest(unittest.makeSuite(Test_Darcy3D))
   suite.addTest(unittest.makeSuite(Test_Darcy2D))
   suite.addTest(unittest.makeSuite(Test_StokesProblemCartesian3D))
   suite.addTest(unittest.makeSuite(Test_Mountains3D))
   suite.addTest(unittest.makeSuite(Test_Mountains2D))
   suite.addTest(unittest.makeSuite(Test_Rheologies))
   suite.addTest(unittest.makeSuite(Test_IncompressibleIsotropicFlowCartesian))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

