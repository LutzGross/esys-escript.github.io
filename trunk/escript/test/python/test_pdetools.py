# $Id$

"""
Test suite for the pdetools module

The tests must be linked with a Domain class object in the setUp method:

   from esys.finley import Rectangle
   class Test_LinearPDEOnFinley(Test_LinearPDE):
       RES_TOL=1-8
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
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

import unittest
from esys.escript import *
from esys.escript.pdetools import Locator,Projector,TimeIntegrationManager,NoPDE

class Test_pdetools(unittest.TestCase):
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

    def testProjector_rank0_fast(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=False,fast=True)
      td_ref=x[0]
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank1_fast(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=False,fast=True)
      td_ref=x
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank2_fast(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=False,fast=True)
      td_ref=[[11.,12.],[21,22.]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank3_fast(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=False,fast=True)
      td_ref=[[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

    def testProjector_rank4_fast(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=False,fast=True)
      td_ref=[[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]]*(x[0]+x[1])
      td=p(td_ref.interpolate(Function(self.domain)))
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

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
      p.setValue(D=x[:2],Y=2*x[:2],q=msk,r=2.)
      u=p.getSolution()
      u_ex=2.
      self.failUnless(Lsup(u_ex-u)<Lsup(u_ex)*self.RES_TOL,"value wrong")
