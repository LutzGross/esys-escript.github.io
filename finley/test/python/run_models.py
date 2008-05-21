#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import unittest
import tempfile
      
from esys.escript import *
from esys.finley import Rectangle
import sys
import os
try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'


NE=6
TOL=1.e-5
VERBOSE=False or True

from esys.escript import *
from esys.escript.models import StokesProblemCartesian
from esys.finley import Rectangle, Brick
class Test_Simple2DModels(unittest.TestCase):
   def setUp(self):
       self.domain=Rectangle(NE,NE,order=2,useFullElementOrder=True)
   def tearDown(self):
       del self.domain
   def test_StokesProblemCartesian_PCG_P_0(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="PCG")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_PCG_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="PCG")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_PCG_P_large(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="PCG")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_GMRES_P_0(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="GMRES",iter_restart=20)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))

      # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_GMRES_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="GMRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_GMRES_P_large(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="GMRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_MINRES_P_0(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="MINRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_MINRES_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="MINRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

#   def test_StokesProblemCartesian_MINRES_P_large(self):
#       ETA=1.
#       P1=1000.
#
#       x=self.domain.getX()
#       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
#       mask=whereZero(x[0])    * [1.,1.] \
#              +whereZero(x[0]-1)  * [1.,1.] \
#              +whereZero(x[1])    * [1.,0.] \
#              +whereZero(x[1]-1)  * [1.,1.]
       
#       sp=StokesProblemCartesian(self.domain)
       
#       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
#       u0=(1-x[0])*x[0]*[0.,1.]
#       p0=Scalar(P1,ReducedSolution(self.domain))
#       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="MINRES")
       
#       error_v0=Lsup(u[0]-u0[0])
#       error_v1=Lsup(u[1]-u0[1])/0.25
#       zz=P1*x[0]*x[1]-p
#       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1
#       self.failUnless(error_p<TOL,"pressure error too large.")
#       self.failUnless(error_v0<TOL,"0-velocity error too large.")
#       self.failUnless(error_v1<TOL,"1-velocity error too large.")


#   def test_StokesProblemCartesian_TFQMR_P_0(self):
#       ETA=1.
#       P1=0.

#       x=self.domain.getX()
#       F=-P1*x[1]*[1.,0]+(2*ETA-P1*x[0])*[0.,1.]
#       mask=whereZero(x[0])    * [1.,1.] \
#              +whereZero(x[0]-1)  * [1.,1.] \
#              +whereZero(x[1])    * [1.,0.] \
#              +whereZero(x[1]-1)  * [1.,1.]
       
#       sp=StokesProblemCartesian(self.domain)
       
#       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
#       u0=(1-x[0])*x[0]*[0.,1.]
#       p0=Scalar(P1,ReducedSolution(self.domain))
#       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="TFQMR")
       
#       error_v0=Lsup(u[0]-u0[0])
#       error_v1=Lsup(u[1]-u0[1])/0.25
#       zz=P1*x[0]*x[1]-p
#       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1
#       self.failUnless(error_p<TOL,"pressure error too large.")
#       self.failUnless(error_v0<TOL,"0-velocity error too large.")
#       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_TFQMR_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="TFQMR")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")

   def test_StokesProblemCartesian_TFQMR_P_large(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="TFQMR")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])/0.25
       zz=P1*x[0]*x[1]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")



class Test_Simple3DModels(unittest.TestCase):
   def setUp(self):
       self.domain=Brick(NE,NE,NE,order=2,useFullElementOrder=True)
   def tearDown(self):
       del self.domain
   def test_StokesProblemCartesian_PCG_P_0(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="PCG")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")

   def test_StokesProblemCartesian_PCG_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="PCG")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")
   def test_StokesProblemCartesian_PCG_P_large(self):
       ETA=1.
       P1=1000.

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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="PCG")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")

   def test_StokesProblemCartesian_GMRES_P_0(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="GMRES",iter_restart=20)
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")
   def test_StokesProblemCartesian_GMRES_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="GMRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")
   def test_StokesProblemCartesian_GMRES_P_large(self):
       ETA=1.
       P1=1000.

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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="GMRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")

   def test_StokesProblemCartesian_MINRES_P_0(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="MINRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")

   def test_StokesProblemCartesian_MINRES_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="MINRES")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")

#   def test_StokesProblemCartesian_MINRES_P_large(self):
#       ETA=1.
#       P1=1000.

#       x=self.domain.getX()
#       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
#       mask=whereZero(x[0])    * [1.,1.,1.] \
#              +whereZero(x[0]-1)  * [1.,1.,1.] \
#              +whereZero(x[1])    * [1.,1.,1.] \
#              +whereZero(x[1]-1)  * [1.,1.,1.] \
#              +whereZero(x[2])    * [1.,1.,0.] \
#              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
#       sp=StokesProblemCartesian(self.domain)
       
#       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
#       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
#       p0=Scalar(P1,ReducedSolution(self.domain))
#       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="MINRES")
       
#       error_v0=Lsup(u[0]-u0[0])
#       error_v1=Lsup(u[1]-u0[1])
#       error_v2=Lsup(u[2]-u0[2])/0.25**2
#       zz=P1*x[0]*x[1]*x[2]-p
#       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
#       self.failUnless(error_p<TOL,"pressure error too large.")
#       self.failUnless(error_v0<TOL,"0-velocity error too large.")
#       self.failUnless(error_v1<TOL,"1-velocity error too large.")
#       self.failUnless(error_v2<TOL,"2-velocity error too large.")

#   def test_StokesProblemCartesian_TFQMR_P_0(self):
#       ETA=1.
#       P1=0.

#       x=self.domain.getX()
#       F=-P1*x[1]*x[2]*[1.,0.,0.]-P1*x[0]*x[2]*[0.,1.,0.]+(2*ETA*((1-x[0])*x[0]+(1-x[1])*x[1])-P1*x[0]*x[1])*[0.,0.,1.]
#       x=self.domain.getX()
#       mask=whereZero(x[0])    * [1.,1.,1.] \
#              +whereZero(x[0]-1)  * [1.,1.,1.] \
#              +whereZero(x[1])    * [1.,1.,1.] \
#              +whereZero(x[1]-1)  * [1.,1.,1.] \
#              +whereZero(x[2])    * [1.,1.,0.] \
#              +whereZero(x[2]-1)  * [1.,1.,1.]
       
       
#       sp=StokesProblemCartesian(self.domain)
       
#       sp.initialize(f=F,fixed_u_mask=mask,eta=ETA)
#       u0=(1-x[0])*x[0]*(1-x[1])*x[1]*[0.,0.,1.]
#       p0=Scalar(P1,ReducedSolution(self.domain))
#       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="TFQMR")
       
#       error_v0=Lsup(u[0]-u0[0])
#       error_v1=Lsup(u[1]-u0[1])
#       error_v2=Lsup(u[2]-u0[2])/0.25**2
#       zz=P1*x[0]*x[1]*x[2]-p
#       error_p=Lsup(zz-integrate(zz))
       # print error_p, error_v0,error_v1,error_v2
#       self.failUnless(error_p<TOL,"pressure error too large.")
#       self.failUnless(error_v0<TOL,"0-velocity error too large.")
#       self.failUnless(error_v1<TOL,"1-velocity error too large.")
#       self.failUnless(error_v2<TOL,"2-velocity error too large.")

   def test_StokesProblemCartesian_TFQMR_P_small(self):
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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="TFQMR")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")

   def test_StokesProblemCartesian_TFQMR_P_large(self):
       ETA=1.
       P1=1000.

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
       u,p=sp.solve(u0,p0,show_details=VERBOSE, verbose=VERBOSE,max_iter=100,solver="TFQMR")
       
       error_v0=Lsup(u[0]-u0[0])
       error_v1=Lsup(u[1]-u0[1])
       error_v2=Lsup(u[2]-u0[2])/0.25**2
       zz=P1*x[0]*x[1]*x[2]-p
       error_p=Lsup(zz-integrate(zz))/P1
       # print error_p, error_v0,error_v1,error_v2
       self.failUnless(error_p<TOL,"pressure error too large.")
       self.failUnless(error_v0<TOL,"0-velocity error too large.")
       self.failUnless(error_v1<TOL,"1-velocity error too large.")
       self.failUnless(error_v2<TOL,"2-velocity error too large.")


if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_Simple2DModels))
   # suite.addTest(unittest.makeSuite(Test_Simple3DModels))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

