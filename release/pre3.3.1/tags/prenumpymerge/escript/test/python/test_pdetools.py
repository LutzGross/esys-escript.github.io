
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

"""
Test suite for the pdetools module

The tests must be linked with a Domain class object in the setUp method:

   from esys.finley import Rectangle
   class Test_LinearPDEOnFinley(Test_LinearPDE):
       RES_TOL=1.e-8
       def setUp(self):
           self.domain = Rectangle(10,10,2)
       def tearDown(self):
           del self.domain
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley))
   unittest.TextTestRunner(verbosity=2).run(suite)

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import unittest
from esys.escript import *
from esys.escript.pdetools import Locator,Projector,TimeIntegrationManager,NoPDE,PCG, ArithmeticTuple, GMRES, MINRES, TFQMR, HomogeneousSaddlePointProblem
from esys.escript.pdetools import Defect, NewtonGMRES
from numarray.linear_algebra import solve_linear_equations

class Test_pdetools_noLumping(unittest.TestCase):
    DEBUG=False
    VERBOSE=False
    def test_TimeIntegrationManager_scalar(self):
        t=0.
        dt=0.1
        tm=TimeIntegrationManager(0.,p=1)
        while t<1.:
           t+=dt
           tm.checkin(dt,t)
        v_guess=tm.extrapolate(dt)
        self.failUnless(abs(v_guess-(tm.getTime()+dt))<self.RES_TOL,"extrapolation is wrong")

    def test_TimeIntegrationManager_vector(self):
        t=0.
        dt=0.3
        tm=TimeIntegrationManager(0.,0.,p=1)
        while t<1.:
           t+=dt
           tm.checkin(dt,t,3*t)
        v_guess=tm.extrapolate(dt)
        e=max(abs(v_guess[0]-(tm.getTime()+dt)),abs(v_guess[1]-(tm.getTime()+dt)*3.))
        self.failUnless(e<self.RES_TOL,"extrapolation is wrong")

    def test_Locator(self):
        x=self.domain.getX()
        l=Locator(self.domain,numarray.ones((self.domain.getDim(),)))
        self.failUnless(ContinuousFunction(self.domain)==l.getFunctionSpace(),"wrong function space from domain")

        l=Locator(ContinuousFunction(self.domain),numarray.ones((self.domain.getDim(),)))
        self.failUnless(ContinuousFunction(self.domain)==l.getFunctionSpace(),"wrong function space")

        xx=l.getX()
        self.failUnless(isinstance(xx,numarray.NumArray),"wrong vector type")
        self.failUnless(Lsup(xx-numarray.ones((self.domain.getDim(),)))<self.RES_TOL,"location wrong")
        xx=l(x)
        self.failUnless(isinstance(xx,numarray.NumArray),"wrong vector type")
        self.failUnless(Lsup(xx-numarray.ones((self.domain.getDim(),)))<self.RES_TOL,"value wrong vector")
        xx=l(x[0]+x[1])
        self.failUnless(isinstance(xx,float),"wrong scalar type")
        self.failUnless(abs(xx-2.)<self.RES_TOL,"value wrong scalar")

    def test_Locator_withList(self):
        x=self.domain.getX()
        arg=[numarray.ones((self.domain.getDim(),)), numarray.zeros((self.domain.getDim(),))]
        l=Locator(self.domain,arg)
        self.failUnless(ContinuousFunction(self.domain)==l.getFunctionSpace(),"wrong function space from domain")

        l=Locator(ContinuousFunction(self.domain),arg)
        self.failUnless(ContinuousFunction(self.domain)==l.getFunctionSpace(),"wrong function space")

        xx=l.getX()
        self.failUnless(isinstance(xx,list),"list expected")
        for i in range(len(xx)):
           self.failUnless(isinstance(xx[i],numarray.NumArray),"vector expected for %s item"%i)
           self.failUnless(Lsup(xx[i]-arg[i])<self.RES_TOL,"%s-th location is wrong"%i)
        xx=l(x)
        self.failUnless(isinstance(xx,list),"list expected (2)")
        for i in range(len(xx)):
           self.failUnless(isinstance(xx[i],numarray.NumArray),"vector expected for %s item (2)"%i)
           self.failUnless(Lsup(xx[i]-arg[i])<self.RES_TOL,"%s-th location is wrong (2)"%i)
        xx=l(x[0]+x[1])
        self.failUnless(isinstance(xx,list),"list expected (3)")
        for i in range(len(xx)):
           self.failUnless(isinstance(xx[i],float),"wrong scalar type")
           self.failUnless(abs(xx[i]-(arg[i][0]+arg[i][1]))<self.RES_TOL,"value wrong scalar")
      
    def testProjector_rank0(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=x[0]
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank1(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=x
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank2(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=[[11.,12.],[21,22.]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank3(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=[[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank4(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=[[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")


    def testProjector_rank0_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=x[0]
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank1_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=x
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank2_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=[[11.,12.],[21,22.]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank3_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=[[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank4_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=[[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank0_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=x[0]
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank1_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=x
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank2_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=[[11.,12.],[21,22.]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank3_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=[[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank4_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=False,fast=False)
      td_ref=[[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")


    def testProjector_rank0_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=1.
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank1_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=numarray.array([1.,2.,3.])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank2_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=numarray.array([[11.,12.],[21,22.]])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank3_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=numarray.array([[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")

    def testProjector_rank4_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      p=Projector(self.domain,reduce=True,fast=False)
      td_ref=numarray.array([[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*self.RES_TOL,"value wrong")


    def test_NoPDE_scalar_missing_r(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])
      p.setValue(D=1.,Y=1.,q=msk)
      u=p.getSolution()
      u_ex=(1.-msk)
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")

    def test_NoPDE_scalar_missing_Y(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])
      p.setValue(D=1.,q=msk,r=2.)
      u=p.getSolution()
      u_ex=msk*2.
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")

    def test_NoPDE_scalar_constant(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])
      p.setValue(D=1.,Y=1.,q=msk,r=2.)
      u=p.getSolution()
      u_ex=(1.-msk)+msk*2.
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")

    def test_NoPDE_scalar_variable(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])
      p.setValue(D=x[0],Y=2*x[0],q=msk,r=2.)
      u=p.getSolution()
      u_ex=2.
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")

    def test_NoPDE_vector_missing_Y(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])*[1.,0.]
      p.setValue(D=numarray.ones([2]),q=msk,r=2.)
      u=p.getSolution()
      u_ex=msk*2.
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")

    def test_NoPDE_vector_missing_r(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])*[1.,0.]
      p.setValue(D=numarray.ones([2]),Y=numarray.ones([2]),q=msk)
      u=p.getSolution()
      u_ex=(1.-msk)
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")

    def test_NoPDE_vector_constant(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])*[1.,0.]
      p.setValue(D=numarray.ones([2]),Y=numarray.ones([2]),q=msk,r=2.)
      u=p.getSolution()
      u_ex=(1.-msk)+msk*2.
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")

    def test_NoPDE_vector_variable(self):
      p=NoPDE(self.domain)
      x=self.domain.getX()
      msk=whereZero(x[0])*[1.,0.]
      p.setValue(D=x[:2]+1,Y=2*(x[:2]+1),q=msk,r=2.)
      u=p.getSolution()
      u_ex=2.
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")
    #=====
    def testPCG(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                  5.834798554135237e-01, -3.704394311709722e+00,
                  5.765369186984777e+00, -1.309786358737351e+01,
                  2.522087134507148e+01, -3.393956279045637e+01,
                  1.046856914770830e+02, -2.447764190849540e+02],
               [ -2.391895572674098e-01,  1.256797283910693e+02,
                 -9.188270412920813e-01,  1.300169538880688e+00,
                 -5.353714719231424e-01,  2.674709444667012e+00,
                 -1.116097841269580e+01,  2.801193427514478e+01,
                 -3.877806125898224e+01,  3.063505753648256e+01],
               [  5.834798554135237e-01, -9.188270412920813e-01,
                  6.240841811806843e+01, -8.176289504109282e-01,
                  1.447935098417076e-01, -9.721424148655324e-01,
                  6.713551574117577e-01, -3.656297654168375e+00,
                  7.015141656913973e+00, -4.195525932156250e+01],
               [ -3.704394311709722e+00,  1.300169538880688e+00,
                 -8.176289504109282e-01,  3.604980536782198e+01,
                 -6.241238423759328e-01,  1.142345320047869e+00,
                 -3.438816797096519e+00,  5.854857481367470e+00,
                 -4.524311288596452e+00,  1.136590280389803e+01],
               [  5.765369186984777e+00, -5.353714719231424e-01,
                  1.447935098417076e-01, -6.241238423759328e-01,
                  2.953997190215862e+01, -9.474729233464712e-01,
                  1.883516378345809e+00, -1.906274765704230e+00,
                  4.401859671778645e+00, -1.064573816075257e+01],
               [ -1.309786358737351e+01,  2.674709444667012e+00,
                 -9.721424148655324e-01,  1.142345320047869e+00,
                 -9.474729233464712e-01,  2.876998216302979e+01,
                 -4.853065259692995e-01,  7.088596468102618e-01,
                 -8.972224295152829e-01,  5.228606946522749e+00],
               [  2.522087134507148e+01, -1.116097841269580e+01,
                  6.713551574117577e-01, -3.438816797096519e+00,
                  1.883516378345809e+00, -4.853065259692995e-01,
                  5.121175860935919e+01, -3.523133115905478e-01,
                  1.782136702229135e+00, -1.560849559916187e+00],
               [ -3.393956279045637e+01,  2.801193427514478e+01,
                 -3.656297654168375e+00,  5.854857481367470e+00,
                 -1.906274765704230e+00,  7.088596468102618e-01,
                 -3.523133115905478e-01,  8.411681423853814e+01,
                 -5.238590858177903e-01,  1.515872114883926e+00],
               [  1.046856914770830e+02, -3.877806125898224e+01,
                  7.015141656913973e+00, -4.524311288596452e+00,
                  4.401859671778645e+00, -8.972224295152829e-01,
                  1.782136702229135e+00, -5.238590858177903e-01,
                  1.797889693808014e+02, -8.362340479938084e-01],
               [ -2.447764190849540e+02,  3.063505753648256e+01,
                 -4.195525932156250e+01,  1.136590280389803e+01,
                 -1.064573816075257e+01,  5.228606946522749e+00,
                 -1.560849559916187e+00,  1.515872114883926e+00,
                 -8.362340479938084e-01,  3.833719335346630e+02]])
      x_ref=array([ 0.41794207085296,   0.031441086046563,  0.882801683420401,
                     0.807186823427233,  0.48950999450145,   0.995486532098031,
                     0.351243009576568,  0.704352576819321,  0.850648989740204,
                     0.314596738052894])
      b=array([ 182.911023960262952,   -1.048322041992754,   44.181293875206201,
                30.344553414038817,   15.247917439094513,   24.060664905403492,
                27.210293789825833,   47.122067744075842,  199.267136417856847,
                -8.7934289814322  ])

      def Ap(x):
          return matrixmultiply(A,x)
      def Ms(b):
          out=zeros((size(b),),Float64)
          for i in xrange(size(b)):
            out[i]=b[i]/A[i,i]
          return out

      tol=1.e-4
      x,r,a_norm=PCG(b*1.,Ap,x_ref*0.,Ms,dot, atol=0, rtol=tol, iter_max=12)
      self.failUnless(Lsup(x-x_ref)<=Lsup(x_ref)*tol*10.,"wrong solution")
      self.failUnless(Lsup(r-(b-matrixmultiply(A,x)))<=Lsup(b)*EPSILON*100.,"wrong solution")

    def testMINRES(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                  5.834798554135237e-01, -3.704394311709722e+00,
                  5.765369186984777e+00, -1.309786358737351e+01,
                  2.522087134507148e+01, -3.393956279045637e+01,
                  1.046856914770830e+02, -2.447764190849540e+02],
               [ -2.391895572674098e-01,  1.256797283910693e+02,
                 -9.188270412920813e-01,  1.300169538880688e+00,
                 -5.353714719231424e-01,  2.674709444667012e+00,
                 -1.116097841269580e+01,  2.801193427514478e+01,
                 -3.877806125898224e+01,  3.063505753648256e+01],
               [  5.834798554135237e-01, -9.188270412920813e-01,
                  6.240841811806843e+01, -8.176289504109282e-01,
                  1.447935098417076e-01, -9.721424148655324e-01,
                  6.713551574117577e-01, -3.656297654168375e+00,
                  7.015141656913973e+00, -4.195525932156250e+01],
               [ -3.704394311709722e+00,  1.300169538880688e+00,
                 -8.176289504109282e-01,  3.604980536782198e+01,
                 -6.241238423759328e-01,  1.142345320047869e+00,
                 -3.438816797096519e+00,  5.854857481367470e+00,
                 -4.524311288596452e+00,  1.136590280389803e+01],
               [  5.765369186984777e+00, -5.353714719231424e-01,
                  1.447935098417076e-01, -6.241238423759328e-01,
                  2.953997190215862e+01, -9.474729233464712e-01,
                  1.883516378345809e+00, -1.906274765704230e+00,
                  4.401859671778645e+00, -1.064573816075257e+01],
               [ -1.309786358737351e+01,  2.674709444667012e+00,
                 -9.721424148655324e-01,  1.142345320047869e+00,
                 -9.474729233464712e-01,  2.876998216302979e+01,
                 -4.853065259692995e-01,  7.088596468102618e-01,
                 -8.972224295152829e-01,  5.228606946522749e+00],
               [  2.522087134507148e+01, -1.116097841269580e+01,
                  6.713551574117577e-01, -3.438816797096519e+00,
                  1.883516378345809e+00, -4.853065259692995e-01,
                  5.121175860935919e+01, -3.523133115905478e-01,
                  1.782136702229135e+00, -1.560849559916187e+00],
               [ -3.393956279045637e+01,  2.801193427514478e+01,
                 -3.656297654168375e+00,  5.854857481367470e+00,
                 -1.906274765704230e+00,  7.088596468102618e-01,
                 -3.523133115905478e-01,  8.411681423853814e+01,
                 -5.238590858177903e-01,  1.515872114883926e+00],
               [  1.046856914770830e+02, -3.877806125898224e+01,
                  7.015141656913973e+00, -4.524311288596452e+00,
                  4.401859671778645e+00, -8.972224295152829e-01,
                  1.782136702229135e+00, -5.238590858177903e-01,
                  1.797889693808014e+02, -8.362340479938084e-01],
               [ -2.447764190849540e+02,  3.063505753648256e+01,
                 -4.195525932156250e+01,  1.136590280389803e+01,
                 -1.064573816075257e+01,  5.228606946522749e+00,
                 -1.560849559916187e+00,  1.515872114883926e+00,
                 -8.362340479938084e-01,  3.833719335346630e+02]])
      x_ref=array([ 0.41794207085296,   0.031441086046563,  0.882801683420401,
                     0.807186823427233,  0.48950999450145,   0.995486532098031,
                     0.351243009576568,  0.704352576819321,  0.850648989740204,
                     0.314596738052894])
      b=array([ 182.911023960262952,   -1.048322041992754,   44.181293875206201,
                30.344553414038817,   15.247917439094513,   24.060664905403492,
                27.210293789825833,   47.122067744075842,  199.267136417856847,
                -8.7934289814322  ])

      def Ap(x):
          return matrixmultiply(A,x)
      def Ms(b):
          out=zeros((size(b),),Float64)
          for i in xrange(size(b)):
            out[i]=b[i]/A[i,i]
          return out

      tol=1.e-4
      x=MINRES(b*1.,Ap,x_ref*0,Ms,dot, atol=0, rtol=tol, iter_max=12)
      self.failUnless(Lsup(x-x_ref)<=Lsup(x_ref)*tol*10.,"wrong solution")

    def testTFQMR(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                  5.834798554135237e-01, -3.704394311709722e+00,
                  5.765369186984777e+00, -1.309786358737351e+01,
                  2.522087134507148e+01, -3.393956279045637e+01,
                  1.046856914770830e+02, -2.447764190849540e+02],
               [ -2.391895572674098e-01,  1.256797283910693e+02,
                 -9.188270412920813e-01,  1.300169538880688e+00,
                 -5.353714719231424e-01,  2.674709444667012e+00,
                 -1.116097841269580e+01,  2.801193427514478e+01,
                 -3.877806125898224e+01,  3.063505753648256e+01],
               [  5.834798554135237e-01, -9.188270412920813e-01,
                  6.240841811806843e+01, -8.176289504109282e-01,
                  1.447935098417076e-01, -9.721424148655324e-01,
                  6.713551574117577e-01, -3.656297654168375e+00,
                  7.015141656913973e+00, -4.195525932156250e+01],
               [ -3.704394311709722e+00,  1.300169538880688e+00,
                 -8.176289504109282e-01,  3.604980536782198e+01,
                 -6.241238423759328e-01,  1.142345320047869e+00,
                 -3.438816797096519e+00,  5.854857481367470e+00,
                 -4.524311288596452e+00,  1.136590280389803e+01],
               [  5.765369186984777e+00, -5.353714719231424e-01,
                  1.447935098417076e-01, -6.241238423759328e-01,
                  2.953997190215862e+01, -9.474729233464712e-01,
                  1.883516378345809e+00, -1.906274765704230e+00,
                  4.401859671778645e+00, -1.064573816075257e+01],
               [ -1.309786358737351e+01,  2.674709444667012e+00,
                 -9.721424148655324e-01,  1.142345320047869e+00,
                 -9.474729233464712e-01,  2.876998216302979e+01,
                 -4.853065259692995e-01,  7.088596468102618e-01,
                 -8.972224295152829e-01,  5.228606946522749e+00],
               [  2.522087134507148e+01, -1.116097841269580e+01,
                  6.713551574117577e-01, -3.438816797096519e+00,
                  1.883516378345809e+00, -4.853065259692995e-01,
                  5.121175860935919e+01, -3.523133115905478e-01,
                  1.782136702229135e+00, -1.560849559916187e+00],
               [ -3.393956279045637e+01,  2.801193427514478e+01,
                 -3.656297654168375e+00,  5.854857481367470e+00,
                 -1.906274765704230e+00,  7.088596468102618e-01,
                 -3.523133115905478e-01,  8.411681423853814e+01,
                 -5.238590858177903e-01,  1.515872114883926e+00],
               [  1.046856914770830e+02, -3.877806125898224e+01,
                  7.015141656913973e+00, -4.524311288596452e+00,
                  4.401859671778645e+00, -8.972224295152829e-01,
                  1.782136702229135e+00, -5.238590858177903e-01,
                  1.797889693808014e+02, -8.362340479938084e-01],
               [ -2.447764190849540e+02,  3.063505753648256e+01,
                 -4.195525932156250e+01,  1.136590280389803e+01,
                 -1.064573816075257e+01,  5.228606946522749e+00,
                 -1.560849559916187e+00,  1.515872114883926e+00,
                 -8.362340479938084e-01,  3.833719335346630e+02]])
      x_ref=array([ 0.41794207085296,   0.031441086046563,  0.882801683420401,
                     0.807186823427233,  0.48950999450145,   0.995486532098031,
                     0.351243009576568,  0.704352576819321,  0.850648989740204,
                     0.314596738052894])
      b=array([ 182.911023960262952,   -1.048322041992754,   44.181293875206201,
                30.344553414038817,   15.247917439094513,   24.060664905403492,
                27.210293789825833,   47.122067744075842,  199.267136417856847,
                -8.7934289814322  ])

      def Ap(x):
          out=matrixmultiply(A,x)
          for i in xrange(size(x)):
            out[i]/=A[i,i]
          return out

      tol=1.e-5
      for i in xrange(size(b)): b[i]/=A[i,i]
      x=TFQMR(b,Ap,x_ref*0,dot, atol=0, rtol=tol, iter_max=12)
      self.failUnless(Lsup(x-x_ref)<=Lsup(x_ref)*tol*10.,"wrong solution")

    def testGMRES(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                  5.834798554135237e-01, -3.704394311709722e+00,
                  5.765369186984777e+00, -1.309786358737351e+01,
                  2.522087134507148e+01, -3.393956279045637e+01,
                  1.046856914770830e+02, -2.447764190849540e+02],
               [ -2.391895572674098e-01,  1.256797283910693e+02,
                 -9.188270412920813e-01,  1.300169538880688e+00,
                 -5.353714719231424e-01,  2.674709444667012e+00,
                 -1.116097841269580e+01,  2.801193427514478e+01,
                 -3.877806125898224e+01,  3.063505753648256e+01],
               [  5.834798554135237e-01, -9.188270412920813e-01,
                  6.240841811806843e+01, -8.176289504109282e-01,
                  1.447935098417076e-01, -9.721424148655324e-01,
                  6.713551574117577e-01, -3.656297654168375e+00,
                  7.015141656913973e+00, -4.195525932156250e+01],
               [ -3.704394311709722e+00,  1.300169538880688e+00,
                 -8.176289504109282e-01,  3.604980536782198e+01,
                 -6.241238423759328e-01,  1.142345320047869e+00,
                 -3.438816797096519e+00,  5.854857481367470e+00,
                 -4.524311288596452e+00,  1.136590280389803e+01],
               [  5.765369186984777e+00, -5.353714719231424e-01,
                  1.447935098417076e-01, -6.241238423759328e-01,
                  2.953997190215862e+01, -9.474729233464712e-01,
                  1.883516378345809e+00, -1.906274765704230e+00,
                  4.401859671778645e+00, -1.064573816075257e+01],
               [ -1.309786358737351e+01,  2.674709444667012e+00,
                 -9.721424148655324e-01,  1.142345320047869e+00,
                 -9.474729233464712e-01,  2.876998216302979e+01,
                 -4.853065259692995e-01,  7.088596468102618e-01,
                 -8.972224295152829e-01,  5.228606946522749e+00],
               [  2.522087134507148e+01, -1.116097841269580e+01,
                  6.713551574117577e-01, -3.438816797096519e+00,
                  1.883516378345809e+00, -4.853065259692995e-01,
                  5.121175860935919e+01, -3.523133115905478e-01,
                  1.782136702229135e+00, -1.560849559916187e+00],
               [ -3.393956279045637e+01,  2.801193427514478e+01,
                 -3.656297654168375e+00,  5.854857481367470e+00,
                 -1.906274765704230e+00,  7.088596468102618e-01,
                 -3.523133115905478e-01,  8.411681423853814e+01,
                 -5.238590858177903e-01,  1.515872114883926e+00],
               [  1.046856914770830e+02, -3.877806125898224e+01,
                  7.015141656913973e+00, -4.524311288596452e+00,
                  4.401859671778645e+00, -8.972224295152829e-01,
                  1.782136702229135e+00, -5.238590858177903e-01,
                  1.797889693808014e+02, -8.362340479938084e-01],
               [ -2.447764190849540e+02,  3.063505753648256e+01,
                 -4.195525932156250e+01,  1.136590280389803e+01,
                 -1.064573816075257e+01,  5.228606946522749e+00,
                 -1.560849559916187e+00,  1.515872114883926e+00,
                 -8.362340479938084e-01,  3.833719335346630e+02]])
      x_ref=array([ 0.41794207085296,   0.031441086046563,  0.882801683420401,
                     0.807186823427233,  0.48950999450145,   0.995486532098031,
                     0.351243009576568,  0.704352576819321,  0.850648989740204,
                     0.314596738052894])
      b=array([ 182.911023960262952,   -1.048322041992754,   44.181293875206201,
                30.344553414038817,   15.247917439094513,   24.060664905403492,
                27.210293789825833,   47.122067744075842,  199.267136417856847,
                -8.7934289814322  ])

      def Ap(x):
          b=matrixmultiply(A,x)
          for i in xrange(size(b)):
            b[i]/=A[i,i]
          return b

      tol=1.e-4
      for i in xrange(size(b)): b[i]/=A[i,i]
      x=GMRES(b,Ap,x_ref*0,dot,atol=0, rtol=tol, iter_max=12)
      self.failUnless(Lsup(x-x_ref)<=Lsup(x_ref)*tol*10.,"wrong solution")

    def testGMRES_P_R(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                  5.834798554135237e-01, -3.704394311709722e+00,
                  5.765369186984777e+00, -1.309786358737351e+01,
                  2.522087134507148e+01, -3.393956279045637e+01,
                  1.046856914770830e+02, -2.447764190849540e+02],
               [ -2.391895572674098e-01,  1.256797283910693e+02,
                 -9.188270412920813e-01,  1.300169538880688e+00,
                 -5.353714719231424e-01,  2.674709444667012e+00,
                 -1.116097841269580e+01,  2.801193427514478e+01,
                 -3.877806125898224e+01,  3.063505753648256e+01],
               [  5.834798554135237e-01, -9.188270412920813e-01,
                  6.240841811806843e+01, -8.176289504109282e-01,
                  1.447935098417076e-01, -9.721424148655324e-01,
                  6.713551574117577e-01, -3.656297654168375e+00,
                  7.015141656913973e+00, -4.195525932156250e+01],
               [ -3.704394311709722e+00,  1.300169538880688e+00,
                 -8.176289504109282e-01,  3.604980536782198e+01,
                 -6.241238423759328e-01,  1.142345320047869e+00,
                 -3.438816797096519e+00,  5.854857481367470e+00,
                 -4.524311288596452e+00,  1.136590280389803e+01],
               [  5.765369186984777e+00, -5.353714719231424e-01,
                  1.447935098417076e-01, -6.241238423759328e-01,
                  2.953997190215862e+01, -9.474729233464712e-01,
                  1.883516378345809e+00, -1.906274765704230e+00,
                  4.401859671778645e+00, -1.064573816075257e+01],
               [ -1.309786358737351e+01,  2.674709444667012e+00,
                 -9.721424148655324e-01,  1.142345320047869e+00,
                 -9.474729233464712e-01,  2.876998216302979e+01,
                 -4.853065259692995e-01,  7.088596468102618e-01,
                 -8.972224295152829e-01,  5.228606946522749e+00],
               [  2.522087134507148e+01, -1.116097841269580e+01,
                  6.713551574117577e-01, -3.438816797096519e+00,
                  1.883516378345809e+00, -4.853065259692995e-01,
                  5.121175860935919e+01, -3.523133115905478e-01,
                  1.782136702229135e+00, -1.560849559916187e+00],
               [ -3.393956279045637e+01,  2.801193427514478e+01,
                 -3.656297654168375e+00,  5.854857481367470e+00,
                 -1.906274765704230e+00,  7.088596468102618e-01,
                 -3.523133115905478e-01,  8.411681423853814e+01,
                 -5.238590858177903e-01,  1.515872114883926e+00],
               [  1.046856914770830e+02, -3.877806125898224e+01,
                  7.015141656913973e+00, -4.524311288596452e+00,
                  4.401859671778645e+00, -8.972224295152829e-01,
                  1.782136702229135e+00, -5.238590858177903e-01,
                  1.797889693808014e+02, -8.362340479938084e-01],
               [ -2.447764190849540e+02,  3.063505753648256e+01,
                 -4.195525932156250e+01,  1.136590280389803e+01,
                 -1.064573816075257e+01,  5.228606946522749e+00,
                 -1.560849559916187e+00,  1.515872114883926e+00,
                 -8.362340479938084e-01,  3.833719335346630e+02]])
      x_ref=array([ 0.41794207085296,   0.031441086046563,  0.882801683420401,
                     0.807186823427233,  0.48950999450145,   0.995486532098031,
                     0.351243009576568,  0.704352576819321,  0.850648989740204,
                     0.314596738052894])
      b=array([ 182.911023960262952,   -1.048322041992754,   44.181293875206201,
                30.344553414038817,   15.247917439094513,   24.060664905403492,
                27.210293789825833,   47.122067744075842,  199.267136417856847,
                -8.7934289814322  ])

      def Ap(x):
          return matrixmultiply(A,x)
      def P_Rp(x):
          out=zeros(size(x), Float64)
          for i in xrange(size(x)):
            out[i]=x[i]/A[i,i]
          return out
       
      tol=1.e-4
      x=GMRES(b,Ap,x_ref*0,dot,atol=0, rtol=tol, iter_max=12,P_R=P_Rp)
      self.failUnless(Lsup(x-x_ref)<=Lsup(x_ref)*tol*10.,"wrong solution")

    def testNewtonGMRES(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      class LL(Defect):
           def __init__(self,*kwargs):
                super(LL, self).__init__(*kwargs)
                self.A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                            5.834798554135237e-01, -3.704394311709722e+00,
                            5.765369186984777e+00, -1.309786358737351e+01,
                            2.522087134507148e+01, -3.393956279045637e+01,
                            1.046856914770830e+02, -2.447764190849540e+02],
                         [ -2.391895572674098e-01,  1.256797283910693e+02,
                           -9.188270412920813e-01,  1.300169538880688e+00,
                           -5.353714719231424e-01,  2.674709444667012e+00,
                           -1.116097841269580e+01,  2.801193427514478e+01,
                           -3.877806125898224e+01,  3.063505753648256e+01],
                         [  5.834798554135237e-01, -9.188270412920813e-01,
                            6.240841811806843e+01, -8.176289504109282e-01,
                            1.447935098417076e-01, -9.721424148655324e-01,
                            6.713551574117577e-01, -3.656297654168375e+00,
                            7.015141656913973e+00, -4.195525932156250e+01],
                         [ -3.704394311709722e+00,  1.300169538880688e+00,
                           -8.176289504109282e-01,  3.604980536782198e+01,
                           -6.241238423759328e-01,  1.142345320047869e+00,
                           -3.438816797096519e+00,  5.854857481367470e+00,
                           -4.524311288596452e+00,  1.136590280389803e+01],
                         [  5.765369186984777e+00, -5.353714719231424e-01,
                            1.447935098417076e-01, -6.241238423759328e-01,
                            2.953997190215862e+01, -9.474729233464712e-01,
                            1.883516378345809e+00, -1.906274765704230e+00,
                            4.401859671778645e+00, -1.064573816075257e+01],
                         [ -1.309786358737351e+01,  2.674709444667012e+00,
                           -9.721424148655324e-01,  1.142345320047869e+00,
                           -9.474729233464712e-01,  2.876998216302979e+01,
                           -4.853065259692995e-01,  7.088596468102618e-01,
                           -8.972224295152829e-01,  5.228606946522749e+00],
                         [  2.522087134507148e+01, -1.116097841269580e+01,
                            6.713551574117577e-01, -3.438816797096519e+00,
                            1.883516378345809e+00, -4.853065259692995e-01,
                            5.121175860935919e+01, -3.523133115905478e-01,
                            1.782136702229135e+00, -1.560849559916187e+00],
                         [ -3.393956279045637e+01,  2.801193427514478e+01,
                           -3.656297654168375e+00,  5.854857481367470e+00,
                           -1.906274765704230e+00,  7.088596468102618e-01,
                           -3.523133115905478e-01,  8.411681423853814e+01,
                           -5.238590858177903e-01,  1.515872114883926e+00],
                         [  1.046856914770830e+02, -3.877806125898224e+01,
                            7.015141656913973e+00, -4.524311288596452e+00,
                            4.401859671778645e+00, -8.972224295152829e-01,
                            1.782136702229135e+00, -5.238590858177903e-01,
                            1.797889693808014e+02, -8.362340479938084e-01],
                         [ -2.447764190849540e+02,  3.063505753648256e+01,
                           -4.195525932156250e+01,  1.136590280389803e+01,
                           -1.064573816075257e+01,  5.228606946522749e+00,
                           -1.560849559916187e+00,  1.515872114883926e+00,
                           -8.362340479938084e-01,  3.833719335346630e+02]])
                self.x_ref=array([ 0.41794207085296,   0.031441086046563,  0.882801683420401,
                     0.807186823427233,  0.48950999450145,   0.995486532098031,
                     0.351243009576568,  0.704352576819321,  0.850648989740204,
                     0.314596738052894])
                self.b=array([ 182.911023960262952,   -1.048322041992754,   44.181293875206201,
                     30.344553414038817,   15.247917439094513,   24.060664905403492,
                     27.210293789825833,   47.122067744075842,  199.267136417856847,
                     -8.7934289814322  ])
           def eval(self,x):
              out=matrixmultiply(self.A,x)-self.b
              for i in xrange(size(self.b)):
                out[i]/=self.A[i,i]
              return out
           def bilinearform(self,x0,x1):
              return dot(x0,x1)
              
      tol=1.e-8
      ll=LL()
      x=NewtonGMRES(LL(),ll.x_ref*0., iter_max=100, sub_iter_max=20, atol=0,rtol=tol, verbose=self.VERBOSE)
      self.failUnless(Lsup(x-ll.x_ref)<=Lsup(ll.x_ref)*tol*10.,"wrong solution")

    def testNewtonGMRES(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      class LL(Defect):
           def __init__(self,*kwargs):
                super(LL, self).__init__(*kwargs)
                self.A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                            5.834798554135237e-01, -3.704394311709722e+00,
                            5.765369186984777e+00, -1.309786358737351e+01,
                            2.522087134507148e+01, -3.393956279045637e+01,
                            1.046856914770830e+02, -2.447764190849540e+02],
                         [ -2.391895572674098e-01,  1.256797283910693e+02,
                           -9.188270412920813e-01,  1.300169538880688e+00,
                           -5.353714719231424e-01,  2.674709444667012e+00,
                           -1.116097841269580e+01,  2.801193427514478e+01,
                           -3.877806125898224e+01,  3.063505753648256e+01],
                         [  5.834798554135237e-01, -9.188270412920813e-01,
                            6.240841811806843e+01, -8.176289504109282e-01,
                            1.447935098417076e-01, -9.721424148655324e-01,
                            6.713551574117577e-01, -3.656297654168375e+00,
                            7.015141656913973e+00, -4.195525932156250e+01],
                         [ -3.704394311709722e+00,  1.300169538880688e+00,
                           -8.176289504109282e-01,  3.604980536782198e+01,
                           -6.241238423759328e-01,  1.142345320047869e+00,
                           -3.438816797096519e+00,  5.854857481367470e+00,
                           -4.524311288596452e+00,  1.136590280389803e+01],
                         [  5.765369186984777e+00, -5.353714719231424e-01,
                            1.447935098417076e-01, -6.241238423759328e-01,
                            2.953997190215862e+01, -9.474729233464712e-01,
                            1.883516378345809e+00, -1.906274765704230e+00,
                            4.401859671778645e+00, -1.064573816075257e+01],
                         [ -1.309786358737351e+01,  2.674709444667012e+00,
                           -9.721424148655324e-01,  1.142345320047869e+00,
                           -9.474729233464712e-01,  2.876998216302979e+01,
                           -4.853065259692995e-01,  7.088596468102618e-01,
                           -8.972224295152829e-01,  5.228606946522749e+00],
                         [  2.522087134507148e+01, -1.116097841269580e+01,
                            6.713551574117577e-01, -3.438816797096519e+00,
                            1.883516378345809e+00, -4.853065259692995e-01,
                            5.121175860935919e+01, -3.523133115905478e-01,
                            1.782136702229135e+00, -1.560849559916187e+00],
                         [ -3.393956279045637e+01,  2.801193427514478e+01,
                           -3.656297654168375e+00,  5.854857481367470e+00,
                           -1.906274765704230e+00,  7.088596468102618e-01,
                           -3.523133115905478e-01,  8.411681423853814e+01,
                           -5.238590858177903e-01,  1.515872114883926e+00],
                         [  1.046856914770830e+02, -3.877806125898224e+01,
                            7.015141656913973e+00, -4.524311288596452e+00,
                            4.401859671778645e+00, -8.972224295152829e-01,
                            1.782136702229135e+00, -5.238590858177903e-01,
                            1.797889693808014e+02, -8.362340479938084e-01],
                         [ -2.447764190849540e+02,  3.063505753648256e+01,
                           -4.195525932156250e+01,  1.136590280389803e+01,
                           -1.064573816075257e+01,  5.228606946522749e+00,
                           -1.560849559916187e+00,  1.515872114883926e+00,
                           -8.362340479938084e-01,  3.833719335346630e+02]])
                self.x_ref=array([ 0.41794207085296,   0.031441086046563,  0.882801683420401,
                     0.807186823427233,  0.48950999450145,   0.995486532098031,
                     0.351243009576568,  0.704352576819321,  0.850648989740204,
                     0.314596738052894])
                self.b=array([ 182.911023960262952,   -1.048322041992754,   44.181293875206201,
                     30.344553414038817,   15.247917439094513,   24.060664905403492,
                     27.210293789825833,   47.122067744075842,  199.267136417856847,
                     -8.7934289814322  ])
           def eval(self,x):
              out=matrixmultiply(self.A,x)-self.b
              for i in xrange(size(self.b)):
                out[i]/=self.A[i,i]
              return out
           def bilinearform(self,x0,x1):
              return dot(x0,x1)
              
      tol=1.e-8
      ll=LL()
      x=NewtonGMRES(LL(),ll.x_ref*0., iter_max=100, sub_iter_max=20, atol=0,rtol=tol, verbose=self.VERBOSE)
      self.failUnless(Lsup(x-ll.x_ref)<=Lsup(ll.x_ref)*tol*10.,"wrong solution")

    def testHomogeneousSaddlePointProblem_PCG(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      class LL(HomogeneousSaddlePointProblem):
           def initialize(self):
               self.A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                                5.834798554135237e-01, -3.704394311709722e+00,
                                5.765369186984777e+00, -1.309786358737351e+01,
                                2.522087134507148e+01, -3.393956279045637e+01,
                                1.046856914770830e+02, -2.447764190849540e+02],
                             [ -2.391895572674098e-01,  1.256797283910693e+02,
                               -9.188270412920813e-01,  1.300169538880688e+00,
                               -5.353714719231424e-01,  2.674709444667012e+00,
                               -1.116097841269580e+01,  2.801193427514478e+01,
                               -3.877806125898224e+01,  3.063505753648256e+01],
                             [  5.834798554135237e-01, -9.188270412920813e-01,
                                6.240841811806843e+01, -8.176289504109282e-01,
                                1.447935098417076e-01, -9.721424148655324e-01,
                                6.713551574117577e-01, -3.656297654168375e+00,
                                7.015141656913973e+00, -4.195525932156250e+01],
                             [ -3.704394311709722e+00,  1.300169538880688e+00,
                               -8.176289504109282e-01,  3.604980536782198e+01,
                               -6.241238423759328e-01,  1.142345320047869e+00,
                               -3.438816797096519e+00,  5.854857481367470e+00,
                               -4.524311288596452e+00,  1.136590280389803e+01],
                             [  5.765369186984777e+00, -5.353714719231424e-01,
                                1.447935098417076e-01, -6.241238423759328e-01,
                                2.953997190215862e+01, -9.474729233464712e-01,
                                1.883516378345809e+00, -1.906274765704230e+00,
                                4.401859671778645e+00, -1.064573816075257e+01],
                             [ -1.309786358737351e+01,  2.674709444667012e+00,
                               -9.721424148655324e-01,  1.142345320047869e+00,
                               -9.474729233464712e-01,  2.876998216302979e+01,
                               -4.853065259692995e-01,  7.088596468102618e-01,
                               -8.972224295152829e-01,  5.228606946522749e+00],
                             [  2.522087134507148e+01, -1.116097841269580e+01,
                                6.713551574117577e-01, -3.438816797096519e+00,
                                1.883516378345809e+00, -4.853065259692995e-01,
                                5.121175860935919e+01, -3.523133115905478e-01,
                                1.782136702229135e+00, -1.560849559916187e+00],
                             [ -3.393956279045637e+01,  2.801193427514478e+01,
                               -3.656297654168375e+00,  5.854857481367470e+00,
                               -1.906274765704230e+00,  7.088596468102618e-01,
                               -3.523133115905478e-01,  8.411681423853814e+01,
                               -5.238590858177903e-01,  1.515872114883926e+00],
                             [  1.046856914770830e+02, -3.877806125898224e+01,
                                7.015141656913973e+00, -4.524311288596452e+00,
                                4.401859671778645e+00, -8.972224295152829e-01,
                                1.782136702229135e+00, -5.238590858177903e-01,
                                1.797889693808014e+02, -8.362340479938084e-01],
                             [ -2.447764190849540e+02,  3.063505753648256e+01,
                               -4.195525932156250e+01,  1.136590280389803e+01,
                               -1.064573816075257e+01,  5.228606946522749e+00,
                               -1.560849559916187e+00,  1.515872114883926e+00,
                               -8.362340479938084e-01,  3.833719335346630e+02]])
               self.x_ref=array([ 0.100225501676291, -0.308862704993209,  0.064097238997721,
                                  0.253012436539738, -0.346223308561905,  0.2425508275422,
                                 -0.194695862196008,  0.09451439391473,   0.302961126826511,
                                 -0.236043777597633] )

               self.Bt=array([[ 0.01627853113636  ,0.06688235764255 , 0.004870689484614],
                             [ 0.062879587145773 ,0.038798770300146, 0.022155850155616],
                             [ 0.09312121957248  ,0.110244632756116, 0.14053347386784 ],
                             [ 0.059000597728388 ,0.090986953740106, 0.035316011834982],
                             [ 0.091209362659698 ,0.13205572801294 , 0.069462874306956],
                             [ 0.077790176986096 ,0.133626423045765, 0.011149969846981],
                             [ 0.01407283482513  ,0.094910926488907, 0.133498532648644],
                             [ 0.025728916673085 ,0.102542818811672, 0.13657268163218 ],
                             [ 0.071254288170748 ,0.071738715618163, 0.078005951991733],
                             [ 0.049463014576779 ,0.103559223780991, 0.003356415647637]])
               self.p_ref = array([ 2.580984952252628 ,4.054090902056985, 0.935138168128546])

               self.b=array([ 123.322775367582238, -51.556206655564573 , 16.220697868056913,
                              6.512480714694167 , -5.727371407390975  , 4.802494840775022,
                             -4.171606044721161 , -1.862366353566293  ,74.850226163257105,
                             -118.602464657076439])

               self.Sinv=array([[ 9313.705360982807179,-5755.536981691270739,  806.289245589733696],
                                 [-5755.536981691271649, 4606.321002756208145,-1630.50619635660928 ],
                                 [  806.289245589733468,-1630.506196356609053, 2145.65035816388945 ]])
           def inner_pBv(self,p,Bv):
              return dot(p,Bv)
           def Bv(self,v):
              return matrixmultiply(transpose(self.Bt),v)
           def inner_p(self,p0,p1):
              return dot(p0,p1)
           def norm_v(self,v):
               return sqrt(dot(v,v))
           def getV(self,p,v):
               out=self.b-matrixmultiply(self.Bt,p)
               return solve_linear_equations(self.A,out)+v*self.getSubProblemTolerance()
           def norm_Bv(self,Bv):
               return sqrt(dot(Bv,Bv))
           def solve_AinvBt(self,p):
               out=matrixmultiply(self.Bt,p)
               return solve_linear_equations(self.A,out)
           def solve_prec(self,Bv):
               out=Bv*1.
               for i in xrange(size(out)): out[i]*=self.Sinv[i,i]
               return out

      tol=1.e-8
      ll=LL()
      ll.initialize()
      ll.setTolerance(tol)
      # ll.setSubToleranceReductionFactor(0.1)
      x,p=ll.solve(ll.x_ref*1.20,ll.p_ref*(-2),max_iter=20, verbose=False, show_details=False, usePCG=True, iter_restart=20,max_correction_steps=3)
      self.failUnless(Lsup(x-ll.x_ref)<=Lsup(ll.x_ref)*tol*10.,"wrong x solution")
      self.failUnless(Lsup(p-ll.p_ref)<=Lsup(ll.p_ref)*tol*10.,"wrong p solution")

    def testHomogeneousSaddlePointProblem_GMRES(self):
      from numarray import array,matrixmultiply, zeros, dot, size, Float64
      from math import sqrt
      class LL(HomogeneousSaddlePointProblem):
           def initialize(self):
               self.A=array([[  4.752141253159452e+02, -2.391895572674098e-01,
                                5.834798554135237e-01, -3.704394311709722e+00,
                                5.765369186984777e+00, -1.309786358737351e+01,
                                2.522087134507148e+01, -3.393956279045637e+01,
                                1.046856914770830e+02, -2.447764190849540e+02],
                             [ -2.391895572674098e-01,  1.256797283910693e+02,
                               -9.188270412920813e-01,  1.300169538880688e+00,
                               -5.353714719231424e-01,  2.674709444667012e+00,
                               -1.116097841269580e+01,  2.801193427514478e+01,
                               -3.877806125898224e+01,  3.063505753648256e+01],
                             [  5.834798554135237e-01, -9.188270412920813e-01,
                                6.240841811806843e+01, -8.176289504109282e-01,
                                1.447935098417076e-01, -9.721424148655324e-01,
                                6.713551574117577e-01, -3.656297654168375e+00,
                                7.015141656913973e+00, -4.195525932156250e+01],
                             [ -3.704394311709722e+00,  1.300169538880688e+00,
                               -8.176289504109282e-01,  3.604980536782198e+01,
                               -6.241238423759328e-01,  1.142345320047869e+00,
                               -3.438816797096519e+00,  5.854857481367470e+00,
                               -4.524311288596452e+00,  1.136590280389803e+01],
                             [  5.765369186984777e+00, -5.353714719231424e-01,
                                1.447935098417076e-01, -6.241238423759328e-01,
                                2.953997190215862e+01, -9.474729233464712e-01,
                                1.883516378345809e+00, -1.906274765704230e+00,
                                4.401859671778645e+00, -1.064573816075257e+01],
                             [ -1.309786358737351e+01,  2.674709444667012e+00,
                               -9.721424148655324e-01,  1.142345320047869e+00,
                               -9.474729233464712e-01,  2.876998216302979e+01,
                               -4.853065259692995e-01,  7.088596468102618e-01,
                               -8.972224295152829e-01,  5.228606946522749e+00],
                             [  2.522087134507148e+01, -1.116097841269580e+01,
                                6.713551574117577e-01, -3.438816797096519e+00,
                                1.883516378345809e+00, -4.853065259692995e-01,
                                5.121175860935919e+01, -3.523133115905478e-01,
                                1.782136702229135e+00, -1.560849559916187e+00],
                             [ -3.393956279045637e+01,  2.801193427514478e+01,
                               -3.656297654168375e+00,  5.854857481367470e+00,
                               -1.906274765704230e+00,  7.088596468102618e-01,
                               -3.523133115905478e-01,  8.411681423853814e+01,
                               -5.238590858177903e-01,  1.515872114883926e+00],
                             [  1.046856914770830e+02, -3.877806125898224e+01,
                                7.015141656913973e+00, -4.524311288596452e+00,
                                4.401859671778645e+00, -8.972224295152829e-01,
                                1.782136702229135e+00, -5.238590858177903e-01,
                                1.797889693808014e+02, -8.362340479938084e-01],
                             [ -2.447764190849540e+02,  3.063505753648256e+01,
                               -4.195525932156250e+01,  1.136590280389803e+01,
                               -1.064573816075257e+01,  5.228606946522749e+00,
                               -1.560849559916187e+00,  1.515872114883926e+00,
                               -8.362340479938084e-01,  3.833719335346630e+02]])
               self.x_ref=array([ 0.100225501676291, -0.308862704993209,  0.064097238997721,
                                  0.253012436539738, -0.346223308561905,  0.2425508275422,
                                 -0.194695862196008,  0.09451439391473,   0.302961126826511,
                                 -0.236043777597633] )

               self.Bt=array([[ 0.01627853113636  ,0.06688235764255 , 0.004870689484614],
                             [ 0.062879587145773 ,0.038798770300146, 0.022155850155616],
                             [ 0.09312121957248  ,0.110244632756116, 0.14053347386784 ],
                             [ 0.059000597728388 ,0.090986953740106, 0.035316011834982],
                             [ 0.091209362659698 ,0.13205572801294 , 0.069462874306956],
                             [ 0.077790176986096 ,0.133626423045765, 0.011149969846981],
                             [ 0.01407283482513  ,0.094910926488907, 0.133498532648644],
                             [ 0.025728916673085 ,0.102542818811672, 0.13657268163218 ],
                             [ 0.071254288170748 ,0.071738715618163, 0.078005951991733],
                             [ 0.049463014576779 ,0.103559223780991, 0.003356415647637]])
               self.p_ref = array([ 2.580984952252628 ,4.054090902056985, 0.935138168128546])

               self.b=array([ 123.322775367582238, -51.556206655564573 , 16.220697868056913,
                              6.512480714694167 , -5.727371407390975  , 4.802494840775022,
                             -4.171606044721161 , -1.862366353566293  ,74.850226163257105,
                             -118.602464657076439])

               self.Sinv=array([[ 9313.705360982807179,-5755.536981691270739,  806.289245589733696],
                                 [-5755.536981691271649, 4606.321002756208145,-1630.50619635660928 ],
                                 [  806.289245589733468,-1630.506196356609053, 2145.65035816388945 ]])
           def inner_pBv(self,p,Bv):
              return dot(p,Bv)
           def Bv(self,v):
              return matrixmultiply(transpose(self.Bt),v)
           def inner_p(self,p0,p1):
              return dot(p0,p1)
           def norm_v(self,v):
               return sqrt(dot(v,v))
           def getV(self,p,v):
               out=self.b-matrixmultiply(self.Bt,p)
               return solve_linear_equations(self.A,out)+v*self.getSubProblemTolerance()
           def norm_Bv(self,Bv):
               return sqrt(dot(Bv,Bv))
           def solve_AinvBt(self,p):
               out=matrixmultiply(self.Bt,p)
               return solve_linear_equations(self.A,out)
           def solve_prec(self,Bv):
               out=Bv*1.
               for i in xrange(size(out)): out[i]*=self.Sinv[i,i]
               return out
              
      tol=1.e-8
      ll=LL()
      ll.initialize()
      ll.setTolerance(tol)
      # ll.setSubToleranceReductionFactor(0.1)
      x,p=ll.solve(ll.x_ref*1.20,ll.p_ref*(-2),max_iter=20, verbose=False, show_details=False, usePCG=False, iter_restart=20,max_correction_steps=3)
      self.failUnless(Lsup(x-ll.x_ref)<=Lsup(ll.x_ref)*tol*10.,"wrong x solution")
      self.failUnless(Lsup(p-ll.p_ref)<=Lsup(ll.p_ref)*tol*10.,"wrong p solution")

    def testArithmeticTuple(self):
        a=ArithmeticTuple(1.,2.)
        self.failUnless(len(a)==2,"wrong length")
        self.failUnless(a[0]==1.,"wrong first item")
        self.failUnless(a[1]==2.,"wrong second item")
        c=a*6.
        self.failUnless(isinstance(c,ArithmeticTuple),"c is not an instance of ArithmeticTuple")
        self.failUnless(len(c)==2,"c has wrong length")
        self.failUnless(c[0]==6.,"c has wrong first item")
        self.failUnless(c[1]==12.,"c has wrong second item")
        b=5.*a
        self.failUnless(isinstance(b,ArithmeticTuple),"b is not an instance of ArithmeticTuple")
        self.failUnless(len(b)==2,"b has wrong length")
        self.failUnless(b[0]==5.,"b has wrong first item")
        self.failUnless(b[1]==10.,"b has wrong second item")
        a+=ArithmeticTuple(3.,4.)
        self.failUnless(a[0]==4.,"wrong first item of inplace update")
        self.failUnless(a[1]==6.,"wrong second item of inplace update")



class Test_pdetools(Test_pdetools_noLumping):
    def testProjector_rank0_fast_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=x[0]
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank1_fast_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=x
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank2_fast_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=[[11.,12.],[21,22.]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank3_fast_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=[[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank4_fast_reduced(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=[[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank0_fast_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=1.
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank1_fast_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=numarray.array([1.,2.,3.])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank2_fast_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=numarray.array([[11.,12.],[21,22.]])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank3_fast_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=numarray.array([[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank4_fast_reduced_with_reduced_input(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=True,fast=True)
      td_ref=numarray.array([[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]])
      td=p(Data(td_ref,ReducedFunction(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

