# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for linearPDEs class

"""

__author__="Jaco du Plessis"

from esys.escript import NonlinearPDE, Symbol, whereZero, grad, sin, cos, symmetric, matrixmult, FunctionOnBoundary
import numpy
import unittest
from esys.escript.linearPDEs import IllegalCoefficient,IllegalCoefficientValue
import numpy as np
from esys.escript.pdetools import Locator


class Test_nonLinearPDEs(unittest.TestCase):
    DEBUG=False
    VERBOSE=False
    
    
    
class Test_nlpde(Test_nonLinearPDEs):
    def test_run(self):
        #test just to confirm nlpde works   
        u=Symbol('u', dim=self.domain.getDim())
        nlpde = NonlinearPDE(self.domain, u)
        x=self.domain.getX()
        gammaD=whereZero(x[0])+whereZero(x[1])
        nlpde.setValue(X=grad(u), Y=5*u, q=gammaD, r=1)
        v=nlpde.getSolution(u=1)
    def test_setVals1eq(self):
        #test setting Coefficients with 1 equation
        dim=self.domain.getDim()
        u=Symbol('u', dim=dim)
        nlpde = NonlinearPDE(self.domain, u)
        x=self.domain.getX()
        gammaD=whereZero(x[0])+whereZero(x[1])
        nlpde.setValue(X=grad(u), Y=5*u, q=gammaD, r=1)
        A=nlpde.getCoefficient("A")
        B=nlpde.getCoefficient("B")
        C=nlpde.getCoefficient("C")
        D=nlpde.getCoefficient("D")
        if dim==2:
            ATest=numpy.empty((2,2), dtype=object)
            ATest[0]=1,0
            ATest[1]=0,1
            BTest=numpy.empty((2,), dtype=object)
            BTest[0]=0
            BTest[1]=0
            CTest=BTest
            self.assertTrue(A-ATest==Symbol(numpy.zeros((2,2))))
            self.assertTrue(B-BTest==Symbol(numpy.zeros((2,))))
            self.assertTrue(C-CTest==Symbol(numpy.zeros((2,))))
            temp=Symbol('temp') 
            self.assertTrue(D-temp.subs(temp,5)==temp.subs(temp,0))
        else:
            ATest=numpy.empty((3,3), dtype=object)
            ATest[0]=1,0,0
            ATest[1]=0,1,0
            ATest[2]=0,0,1
            BTest=numpy.empty((3,), dtype=object)
            BTest[0]=0
            BTest[1]=0
            BTest[2]=0
            CTest=BTest
            self.assertTrue(A-ATest==Symbol(numpy.zeros((3,3))))
            self.assertTrue(B-BTest==Symbol(numpy.zeros((3,))))
            self.assertTrue(C-CTest==Symbol(numpy.zeros((3,))))
            temp=Symbol('temp') 
            self.assertTrue(D-temp.subs(temp,5)==temp.subs(temp,0))

    def test_setVals2eq(self):
        #test setting Coefficients with 2 coeficients
        dim=self.domain.getDim()
        u = Symbol('u',(2,), dim=dim)
        q = Symbol('q', (2,2))
        nlpde = NonlinearPDE(self.domain, u, debug=0)
        x = self.domain.getX()
        gammaD=whereZero(x[1])*[1,1]#+whereZero(x[0])*[1,0]+whereZero(x[0]-1)*[1,0]  
        yconstraint = FunctionOnBoundary(self.domain).getX()[1]
        nlpde.setValue(X=grad(u),q=gammaD,Y=[-50,-50]*u)
        A=nlpde.getCoefficient("A")
        B=nlpde.getCoefficient("B")
        C=nlpde.getCoefficient("C")
        D=nlpde.getCoefficient("D")
        if dim==2:
            ATest=numpy.empty((2,2,2,2),dtype=object)
            ATest[0]=(((1,0),(0,0)),((0,1),(0,0)))
            ATest[1]=(((0,0),(1,0)),((0,0),(0,1)))
            BTest=numpy.empty((2,2,2),dtype=object)
            BTest[0]=0
            BTest[1]=0
            CTest=BTest
            DTest=numpy.empty((2,2),dtype=object)
            DTest[0]=(-50,0)
            DTest[1]=(0,-50)
            self.assertTrue(numpy.ndarray.__eq__(ATest, A).all())
            self.assertTrue(numpy.ndarray.__eq__(BTest, B).all())
            self.assertTrue(numpy.ndarray.__eq__(CTest, C).all())
            self.assertTrue(numpy.ndarray.__eq__(DTest, D).all())
        else:
            ATest=numpy.empty((2,3,2,3),dtype=object)
            ATest[0]=(((1,0,0),(0,0,0)),((0,1,0),(0,0,0)),((0,0,1),(0,0,0)))
            ATest[1]=(((0,0,0),(1,0,0)),((0,0,0),(0,1,0)),((0,0,0),(0,0,1)))
            #ATest[1]=(((0,0,1),(1,0,0)),((0,0),(0,1)))
            BTest=numpy.empty((2,3,2),dtype=object)
            BTest[0]=0
            BTest[1]=0
            CTest=numpy.empty((2,2,3),dtype=object)
            CTest[0]=0
            CTest[1]=0
            DTest=numpy.empty((2,2),dtype=object)
            DTest[0]=(-50,0)
            DTest[1]=(0,-50)
            self.assertTrue(numpy.ndarray.__eq__(BTest, B).all())
            self.assertTrue(numpy.ndarray.__eq__(CTest, C).all())
            self.assertTrue(numpy.ndarray.__eq__(DTest, D).all())

    def test_DimAndShape1eq(self):
        dim=self.domain.getDim()
        if dim==3:
            u = Symbol('u', dim=2)
        else:
            u = Symbol('u', dim=3)
        nlpde = NonlinearPDE(self.domain, u, debug=0)
        args=dict(X=grad(u), Y=5*u)
        self.assertRaises(IllegalCoefficientValue, nlpde.setValue,**args)
        u = Symbol('u', dim=dim)
        args=dict(X=u, Y=5*u)
        self.assertRaises(IllegalCoefficientValue, nlpde.setValue,**args)
        args=dict(X=grad(u), Y=5*grad(u))
        self.assertRaises(IllegalCoefficientValue, nlpde.setValue,**args)
        #args=dict(q=u)
        #self.assertRaises(IllegalCoefficientValue, nlpde.setValue,**args)
    
    def test_DimAndShape2eq(self):
        dim=self.domain.getDim()
        u = Symbol('u',(2,), dim=dim)
        nlpde = NonlinearPDE(self.domain, u, debug=0)
        args=dict(X=grad(u), Y=5*u[0])
        self.assertRaises(IllegalCoefficientValue, nlpde.setValue,**args)
        args=dict(X=grad(u[0]), Y=5*u)
        self.assertRaises(IllegalCoefficientValue, nlpde.setValue,**args)


    def test_setUnknownPeram(self):
        dim=self.domain.getDim()
        u = Symbol('u',(2,), dim=dim)
        nlpde = NonlinearPDE(self.domain, u, debug=0)
        args=dict(k=0,f=8)  
        self.assertRaises(IllegalCoefficient,nlpde.setValue,**args)


    def test_yDirection(self):
        dim=self.domain.getDim()
        u = Symbol('u',(2,), dim=dim)
        q = Symbol('q', (2,2))
        theta = Symbol('theta')
        theta=3.141/6
        q[0,0]=cos(theta)
        q[0,1]=-sin(theta)
        q[1,0]=sin(theta)
        q[1,1]=cos(theta)
        sigma = Symbol('sigma',(2,2))
        p = NonlinearPDE(self.domain, u, debug=0)
        epsilon = symmetric(grad(u))
        #epsilon = matrixmult(matrixmult(q,epsilon0),q.transpose(1))
        c00=10;c01=8;c05=0
        c01=8;c11=10;c15=0
        c05=0;c15=0;c55=1
        sigma[0,0]=c00*epsilon[0,0]+c01*epsilon[1,1]+c05*2*epsilon[1,0]
        sigma[1,1]=c01*epsilon[0,0]+c11*epsilon[1,1]+c15*2*epsilon[1,0]
        sigma[0,1]=c05*epsilon[0,0]+c15*epsilon[1,1]+c55*2*epsilon[1,0]
        sigma[1,0]=sigma[0,1]
        #sigma0=matrixmult(matrixmult(q.transpose(1),epsilon),q)
        x = self.domain.getX()
        gammaD=whereZero(x[1])*[1,1]#+whereZero(x[0])*[1,0]+whereZero(x[0]-1)*[1,0]  
        yconstraint = FunctionOnBoundary(self.domain).getX()[1]
        p.setValue(X=sigma,q=gammaD,y=[-50,0]*whereZero(yconstraint-1),r=[1,1])
        v = p.getSolution(u=[0,0])
        x=np.ndarray((2,))
        x[0]=0.5
        x[1]=0.5
        loc=Locator(v.getFunctionSpace(),x)
        valAtX=loc(v)
        self.assertTrue(valAtX[0]>10*valAtX[1])
