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
