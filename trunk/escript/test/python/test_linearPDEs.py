# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for linearPDEs class

"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from esys.escript.util import Lsup,kronecker,interpolate,whereZero, outer, swap_axes
from esys.escript import Function,FunctionOnBoundary,FunctionOnContactZero,Solution,ReducedSolution,Vector,ContinuousFunction,Scalar, ReducedFunction,ReducedFunctionOnBoundary,ReducedFunctionOnContactZero,Data, Tensor4, Tensor, getEscriptParamInt
from esys.escript.linearPDEs import LinearPDE,IllegalCoefficientValue,Poisson, IllegalCoefficientFunctionSpace, TransportPDE, IllegalCoefficient, Helmholtz, LameEquation, SolverOptions
import numpy
import unittest

class Test_linearPDEs(unittest.TestCase):
    TOL=1.e-6
    SOLVER_TOL=1.e-10
    DEBUG=False
    VERBOSE=False
    _domainCanInterpolateAdvanced=None

    # Can the domain interpolate from ReducedFunction to Function?
    def specialInterpolationSupported(self):
        if self._domainCanInterpolateAdvanced is None:
            d0=Data(0,(),ReducedFunction(self.domain))
            try:
                d1=Data(d0, Function(self.domain))
                self._domainCanInterpolateAdvanced=True
            except e:
                self._domainCanInterpolateAdvanced=False
        return self._domainCanInterpolateAdvanced

    def check(self,arg,ref_arg,tol=None):
        """
        checks if arg and ref_arg are nearly identical using the `Lsup`
        """
        if tol==None: tol=self.TOL
        return Lsup(arg-ref_arg)<=tol*Lsup(ref_arg)
    
class Test_LameEquation(Test_linearPDEs):

    def test_config(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        d=self.domain.getDim()
        self.assertEqual((mypde.getNumEquations(), mypde.getNumSolutions(), mypde.getSolverOptions().isSymmetric()),(d,d,True),"set up incorrect")

    def test_setCoefficient_q(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(q=x)

        q_ref=interpolate(x,Solution(self.domain))
        self.assertTrue(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("q"),q_ref),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")


    def test_setCoefficient_r(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(r=x)

        r_ref=interpolate(x,Solution(self.domain))
        self.assertTrue(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("r"),r_ref),"r is nor x")
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")


    def test_setCoefficient_F(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(F=x)

        Y_ref=interpolate(x,Function(self.domain))
        self.assertTrue(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("Y"),Y_ref),"Y is not x")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_f(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(f=x)

        y_ref=interpolate(x,FunctionOnBoundary(self.domain))
        self.assertTrue(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("y"),y_ref),"d is not x[0]")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")       
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_sigma(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(sigma=outer(x,x))

        X_ref=interpolate(outer(x,x),Function(self.domain))
        self.assertTrue(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("X"),X_ref),"X is not x X x")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_lambda(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(lame_lambda=x[0])


        k3=kronecker(Function(self.domain))
        k3Xk3=outer(k3,k3)
        A_ref=x[0]*k3Xk3

        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")       
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_mu(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(lame_mu=x[0])


        k3=kronecker(Function(self.domain))
        k3Xk3=outer(k3,k3)
        A_ref=x[0]*(swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3))

        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")        
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_lambdamu(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(lame_lambda=x[0], lame_mu=x[1])

        k3=kronecker(Function(self.domain))
        k3Xk3=outer(k3,k3)
        A_ref=x[0]*k3Xk3+x[1]*(swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3))

        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_solve(self):
       d=self.domain.getDim()
       mypde=LameEquation(self.domain)
       cf=ContinuousFunction(self.domain)
       x=cf.getX()
       u_ex=x
       msk=Vector(0.,cf)
       for i in range(d): msk[i]=whereZero(x[i])
       mypde.setValue(q=msk,r=u_ex,lame_mu=3,lame_lambda=50,f=(2*3+50*d)*FunctionOnBoundary(self.domain).getNormal())

       u=mypde.getSolution()
       self.assertTrue(self.check(u,u_ex,10*self.TOL),"incorrect solution")

class Test_Helmholtz(Test_linearPDEs):

    def test_config(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        self.assertEqual((mypde.getNumEquations(), mypde.getNumSolutions(), mypde.getSolverOptions().isSymmetric()),(1,1,True),"set up incorrect")
    def test_setCoefficient_q(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(q=whereZero(x[0]))

        q_ref=interpolate(whereZero(x[0]),Solution(self.domain))
        A_ref=kronecker(self.domain)

        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")       
        self.assertTrue(self.check(mypde.getCoefficient("q"),q_ref),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_r(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(r=x[0])

        r_ref=interpolate(x[0],Solution(self.domain))
        A_ref=kronecker(self.domain)
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty") 
        self.assertTrue(self.check(mypde.getCoefficient("r"),r_ref),"r is nor x[0]")
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")


    def test_setCoefficient_f(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(f=x[0])

        Y_ref=interpolate(x[0],Function(self.domain))
        A_ref=kronecker(self.domain)
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("Y"),Y_ref),"Y is not x[0]")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")       
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_alpha(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(alpha=x[0])

        d_ref=interpolate(x[0],FunctionOnBoundary(self.domain))
        A_ref=kronecker(self.domain)
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("d"),d_ref),"d is not x[0]")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")       
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_g(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(g=x[0])

        y_ref=interpolate(x[0],FunctionOnBoundary(self.domain))
        A_ref=kronecker(self.domain)
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("y"),y_ref),"y is not x[0]")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")       
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_omega(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(omega=x[0])

        D_ref=interpolate(x[0],Function(self.domain))
        A_ref=kronecker(self.domain)
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("D"),D_ref),"D is not x[0]")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")        
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_solve(self):
       d=self.domain.getDim()
       cf=ContinuousFunction(self.domain)
       u_ex=Scalar(1.,cf)
       mypde=Helmholtz(self.domain)
       mypde.setValue(f=3,omega=3,alpha=2,g=2)
       u=mypde.getSolution()
       self.assertTrue(self.check(u,u_ex,10*self.TOL),"incorrect solution")

class Test_Poisson(Test_linearPDEs):

    def test_config(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        self.assertEqual((mypde.getNumEquations(), mypde.getNumSolutions(), mypde.getSolverOptions().isSymmetric()),(1,1,True),"set up incorrect")
    def test_setCoefficient_q(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        q_ref=interpolate(whereZero(x[0]),Solution(self.domain))
        A_ref=kronecker(self.domain)
        mypde.setValue(q=whereZero(x[0]))
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("q"),q_ref),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")
    def test_setCoefficient_f(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        Y_ref=interpolate(x[0],Function(self.domain))
        A_ref=kronecker(self.domain)
        mypde.setValue(f=x[0])
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("Y"),Y_ref),"Y is not x[0]")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")
    def test_setCoefficient_f_reduced(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        Y_ref=interpolate(x[0],ReducedFunction(self.domain))
        A_ref=kronecker(self.domain)
        mypde.setValue(f_reduced=x[0])
        self.assertTrue(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.assertTrue(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.assertTrue(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.assertTrue(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.assertTrue(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.assertTrue(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.assertTrue(self.check(mypde.getCoefficient("Y_reduced"),Y_ref),"Y_reduced is not x[0]")
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.assertTrue(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        if self.domain.supportsContactElements():
            self.assertTrue(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
            self.assertTrue(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")        
        self.assertTrue(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty")
    def test_solve(self):
       d=self.domain.getDim()
       cf=ContinuousFunction(self.domain)
       x=cf.getX()
       #construct exact solution:
       u_ex=Scalar(1.,cf)
       for i in range(d):
         u_ex*=x[i]*(2.-x[i])
       #construct mask:
       msk=Scalar(0.,cf)
       for i in range(d):
         msk+=whereZero(x[i])
       #construct right hand side
       f=Scalar(0,cf)
       for i in range(d):
          f_p=Scalar(1,cf)
          for j in range(d):
             if i==j:
                f_p*=2.
             else:
                f_p*=x[j]*(2-x[j])
          f+=f_p
       mypde=Poisson(self.domain)
       mypde.setValue(f=f,q=msk)
       u=mypde.getSolution()
       self.assertTrue(self.check(u,u_ex,10*self.TOL),"incorrect solution")

class Test_LinearPDE_noLumping(Test_linearPDEs):
    N=4

    def test_SolverOptions(self):
        so=SolverOptions()
        
        self.assertTrue(so.getSmoother() == 28, "initial Smoother is wrong.")
        self.assertRaises(ValueError,so.setSmoother,-1)
        so.setSmoother(so.GAUSS_SEIDEL)
        self.assertTrue(so.getSmoother() == 28, "Gauss-Seidel smoother is not set.")
        so.setSmoother(so.JACOBI)
        self.assertTrue(so.getSmoother() == 10, "Jacobi smoother is not set.")

        self.assertTrue(so.getLevelMax() == 100, "initial  LevelMax is wrong.")
        self.assertRaises(ValueError,so.setLevelMax,-1)
        so.setLevelMax(20)
        self.assertTrue(so.getLevelMax() == 20, "LevelMax is wrong.")

        self.assertTrue(so.getCoarseningThreshold() == 0.25, "initial  CoarseningThreshold is wrong.")
        self.assertRaises(ValueError,so.setCoarseningThreshold,-1)
        so.setCoarseningThreshold(0.1)
        self.assertTrue(so.getCoarseningThreshold() == 0.1, "CoarseningThreshold is wrong.")
        
        self.assertTrue(so.getMinCoarseMatrixSize() == 500, "initial  Minimum Coarse Matrix Size is wrong.")
        self.assertRaises(ValueError,so.setMinCoarseMatrixSize,-1)
        so.setMinCoarseMatrixSize(1000)
        self.assertTrue(so.getMinCoarseMatrixSize() == 1000, "Minimum Coarse Matrix Size is wrong.")

        self.assertTrue(so.getNumSweeps() == 1, "initial  Sweeps is wrong.")
        self.assertRaises(ValueError,so.setNumSweeps,-1)
        so.setNumSweeps(3)
        self.assertTrue(so.getNumSweeps() == 3, "Sweeps is wrong.")

        self.assertTrue(so.getNumPreSweeps() == 1, "initial  PreSweeps is wrong.")
        self.assertRaises(ValueError,so.setNumPreSweeps,-1)
        so.setNumPreSweeps(4)
        self.assertTrue(so.getNumPreSweeps() == 4, "PreSweeps is wrong.")

        self.assertTrue(so.getNumPostSweeps() == 1, "initial  PostSweeps is wrong.")
        self.assertRaises(ValueError,so.setNumPostSweeps,-1)
        so.setNumPostSweeps(5)
        self.assertTrue(so.getNumPostSweeps() == 5, "PostSweeps is wrong.")

        self.assertTrue(so.getTolerance() == 1.e-8, "initial Tolerance is wrong.")
        self.assertRaises(ValueError,so.setTolerance,-1)
        so.setTolerance(0.2)
        self.assertTrue(so.getTolerance() == 0.2, "Tolerance is wrong.")

        self.assertTrue(so.getAbsoluteTolerance() == 0., "initial  AbsoluteTolerance is wrong.")
        self.assertRaises(ValueError,so.setAbsoluteTolerance,-1)
        so.setAbsoluteTolerance(0.3)
        self.assertTrue(so.getAbsoluteTolerance() == 0.3, "AbsoluteTolerance is wrong.")

        self.assertTrue(so.getInnerTolerance() == 0.9, "initial InnerTolerance is wrong.")
        self.assertRaises(ValueError,so.setInnerTolerance,-1)
        so.setInnerTolerance(0.4)
        self.assertTrue(so.getInnerTolerance() == 0.4, "InnerTolerance is wrong.")

        self.assertTrue(so.getDropTolerance() == 0.01, "initial DropTolerance is wrong.")
        self.assertRaises(ValueError,so.setDropTolerance,-1)
        so.setDropTolerance(0.5)
        self.assertTrue(so.getDropTolerance() == 0.5, "DropDropTolerance is wrong.")

        self.assertTrue(so.getDropStorage() == 2., "initial DropStorage is wrong.")
        self.assertRaises(ValueError,so.setDropStorage,-1)
        so.setDropStorage(10)
        self.assertTrue(so.getDropStorage() == 10, "DropStorage is wrong.")
        
        self.assertTrue(so.getRelaxationFactor() == 0.3, "initial RelaxationFactor is wrong.")
        self.assertRaises(ValueError,so.setRelaxationFactor,-1)
        so.setRelaxationFactor(0.1)
        self.assertTrue(so.getRelaxationFactor() == 0.1, "Relaxation is wrong.")


        self.assertTrue(so.getIterMax() == 100000, "initial IterMax is wrong.")
        self.assertRaises(ValueError,so.setIterMax,0)
        so.setIterMax(11)
        self.assertTrue(so.getIterMax() == 11, "IterMax is wrong.")

        self.assertTrue(so.getInnerIterMax() == 10, "initial InnerIterMax is wrong.")
        self.assertRaises(ValueError,so.setInnerIterMax,0)
        so.setInnerIterMax(12)
        self.assertTrue(so.getInnerIterMax() == 12, "InnerIterMax is wrong.")

        self.assertTrue(so.getTruncation() == 20, "initial Truncation is wrong.")
        self.assertRaises(ValueError,so.setTruncation,0)
        so.setTruncation(13)
        self.assertTrue(so.getTruncation() == 13, "Truncation is wrong.")

        self.assertTrue(so.getRestart() == None, "initial Truncation is wrong.")
        self.assertRaises(ValueError,so.setTruncation,0)
        so.setRestart(14)
        self.assertTrue(so.getRestart() == 14, "Truncation is wrong.")
        so.setRestart(None)
        self.assertTrue(so.getRestart() == None, "Truncation is wrong.")
        
        self.assertTrue(not so.isVerbose(), "initial verbosity flag is wrong.")
        so.setVerbosityOn()
        self.assertTrue(so.isVerbose(), "verbosity (1) flag is wrong.")
        so.setVerbosityOff()
        self.assertTrue(not so.isVerbose(), "verbosity (2) flag is wrong.")
        so.setVerbosity(verbose=True)
        self.assertTrue(so.isVerbose(), "verbosity (3) flag is wrong.")
        so.setVerbosity(verbose=False)
        self.assertTrue(not so.isVerbose(), "verbosity (4) flag is wrong.")

        self.assertTrue(not so.isSymmetric(), "initial symmetry flag is wrong.")
        so.setSymmetryOn()
        self.assertTrue(so.isSymmetric(), "symmetry (1) flag is wrong.")
        so.setSymmetryOff()
        self.assertTrue(not so.isSymmetric(), "symmetry (2) flag is wrong.")
        so.setSymmetry(flag=True)
        self.assertTrue(so.isSymmetric(), "symmetry (3) flag is wrong.")
        so.setSymmetry(flag=False)
        self.assertTrue(not so.isSymmetric(), "symmetry (4) flag is wrong.")

        self.assertTrue(so.adaptInnerTolerance(), "initial InnerToleranceAdaption flag is wrong.")
        so.setInnerToleranceAdaptionOn()
        self.assertTrue(so.adaptInnerTolerance(), "InnerToleranceAdaption (1) flag is wrong.")
        so.setInnerToleranceAdaptionOff()
        self.assertTrue(not so.adaptInnerTolerance(), "InnerToleranceAdaption (2) flag is wrong.")
        so.setInnerToleranceAdaption(adapt=True)
        self.assertTrue(so.adaptInnerTolerance(), "InnerToleranceAdaption (3) flag is wrong.")
        so.setInnerToleranceAdaption(adapt=False)
        self.assertTrue(not so.adaptInnerTolerance(), "InnerToleranceAdaption (4) flag is wrong.")
     
        self.assertTrue(not so.acceptConvergenceFailure(), "initial acceptConvergenceFailure flag is wrong.")
        so.setAcceptanceConvergenceFailureOn()
        self.assertTrue(so.acceptConvergenceFailure(), "acceptConvergenceFailure (1) flag is wrong.")
        so.setAcceptanceConvergenceFailureOff()
        self.assertTrue(not so.acceptConvergenceFailure(), "acceptConvergenceFailure (2) flag is wrong.")
        so.setAcceptanceConvergenceFailure(accept=True)
        self.assertTrue(so.acceptConvergenceFailure(), "acceptConvergenceFailure (3) flag is wrong.")
        so.setAcceptanceConvergenceFailure(accept=False)
        self.assertTrue(not so.acceptConvergenceFailure(), "acceptConvergenceFailure (4) flag is wrong.")   
        
        self.assertTrue(so.getReordering() == 30, "initial Reordering is wrong.")
        self.assertRaises(ValueError,so.setReordering,-1)
        so.setReordering(so.NO_REORDERING)
        self.assertTrue(so.getReordering() == 17, "NO_REORDERING is not set.")
        so.setReordering(so.MINIMUM_FILL_IN)
        self.assertTrue(so.getReordering() == 18, "MINIMUM_FILL_IN is not set.")
        so.setReordering(so.NESTED_DISSECTION)
        self.assertTrue(so.getReordering() == 19, "NESTED_DISSECTION is not set.")
        so.setReordering(so.DEFAULT_REORDERING)
        self.assertTrue(so.getReordering() == 30, "DEFAULT_REORDERING is not set.")
        
        self.assertTrue(so.getPackage() == 0, "initial solver package is wrong.")
        self.assertRaises(ValueError,so.setPackage,-1)
        so.setPackage(so.PASO)
        self.assertTrue(so.getPackage() == 21, "PASO is not set.")
        so.setPackage(so.SUPER_LU)
        self.assertTrue(so.getPackage() == 31, "SUPER_LU is not set.")
        so.setPackage(so.PASTIX)
        self.assertTrue(so.getPackage() == 32, "PASTIX is not set.")
        so.setPackage(so.MKL)
        self.assertTrue(so.getPackage() == 15, "MKL is not set.")
        so.setPackage(so.UMFPACK)
        self.assertTrue(so.getPackage() == 16, "UMFPACK is not set.")
        so.setPackage(so.TRILINOS)
        self.assertTrue(so.getPackage() == 24, "TRILINOS is not set.")

        self.assertTrue(so.getSolverMethod() == 0, "initial SolverMethod is wrong.")
        self.assertRaises(ValueError,so.setSolverMethod,-1)
        so.setSolverMethod(so.DIRECT)
        self.assertTrue(so.getSolverMethod() == 1, "DIRECT is not set.")
        so.setSolverMethod(so.CHOLEVSKY)
        self.assertTrue(so.getSolverMethod() == 2, "CHOLEVSKY is not set.")
        so.setSolverMethod(so.PCG)
        self.assertTrue(so.getSolverMethod() == 3, "PCG is not set.")
        so.setSolverMethod(so.CR)
        self.assertTrue(so.getSolverMethod() == 4, "CR is not set.")
        so.setSolverMethod(so.CGS)
        self.assertTrue(so.getSolverMethod() == 5, "CGS is not set.")
        so.setSolverMethod(so.BICGSTAB)
        self.assertTrue(so.getSolverMethod() == 6, "BICGSTAB is not set.")
        so.setSolverMethod(so.GMRES)
        self.assertTrue(so.getSolverMethod() == 11, "GMRES is not set.")
        so.setSolverMethod(so.PRES20)
        self.assertTrue(so.getSolverMethod() == 12, "PRES20 is not set.")
        so.setSolverMethod(so.LUMPING)
        self.assertTrue(so.getSolverMethod() == 13, "LUMPING is not set.")
        so.setSolverMethod(so.ITERATIVE)
        self.assertTrue(so.getSolverMethod() == 20, "ITERATIVE is not set.")
        so.setSolverMethod(so.NONLINEAR_GMRES)
        self.assertTrue(so.getSolverMethod() == 25, "NONLINEAR_GMRES is not set.")
        so.setSolverMethod(so.TFQMR)
        self.assertTrue(so.getSolverMethod() == 26, "TFQMR is not set.")
        so.setSolverMethod(so.MINRES)
        self.assertTrue(so.getSolverMethod() == 27, "MINRES is not set.")
        so.setSolverMethod(so.DEFAULT)
        self.assertTrue(so.getSolverMethod() == 0, "DEFAULT is not set.")

        self.assertTrue(so.getPreconditioner() == 10, "initial Preconditioner is wrong.")
        self.assertRaises(ValueError,so.setPreconditioner,-1)
        so.setPreconditioner(so.ILU0)
        self.assertTrue(so.getPreconditioner() == 8, "ILU0 is not set.")
        so.setPreconditioner(so.ILUT)
        self.assertTrue(so.getPreconditioner() == 9, "ILUT is not set.")
        so.setPreconditioner(so.JACOBI)
        self.assertTrue(so.getPreconditioner() == 10, "JACOBI is not set.")
        if getEscriptParamInt('DISABLE_AMG', 0):
            print("AMG test disabled on MPI build")
        else:
            so.setPreconditioner(so.AMG)
            self.assertTrue(so.getPreconditioner() == 22, "AMG is not set.")
        so.setPreconditioner(so.REC_ILU)
        self.assertTrue(so.getPreconditioner() == 23, "REC_ILU is not set.")
        so.setPreconditioner(so.GAUSS_SEIDEL)
        self.assertTrue(so.getPreconditioner() == 28, "GAUSS_SEIDEL is not set.")
        so.setPreconditioner(so.RILU)
        self.assertTrue(so.getPreconditioner() == 29, "RILU is not set.")
        so.setPreconditioner(so.AMLI)
        self.assertTrue(so.getPreconditioner() == 38, "AMLI is not set.")
        so.setPreconditioner(so.NO_PRECONDITIONER)
        self.assertTrue(so.getPreconditioner() == 36, "NO_PRECONDITIONER is not set.")        

        self.assertTrue(so.getCoarsening() == 0, "initial Coarseningr is wrong.")
        self.assertRaises(ValueError,so.setCoarsening,-1)
        so.setCoarsening(so.YAIR_SHAPIRA_COARSENING)
        self.assertTrue(so.getCoarsening() == 33, "YAIR_SHAPIRA_COARSENING is not set.")
        so.setCoarsening(so.RUGE_STUEBEN_COARSENING)
        self.assertTrue(so.getCoarsening() == 34, "RUGE_STUEBEN_COARSENING is not set.")
        so.setCoarsening(so.AGGREGATION_COARSENING)
        self.assertTrue(so.getCoarsening() == 35, "AGREGATION_COARSENING is not set.")
        so.setCoarsening(so.STANDARD_COARSENING)
        self.assertTrue(so.getCoarsening() == 39, "STANDARD_COARSENING is not set.")
        so.setCoarsening(so.DEFAULT)
        self.assertTrue(so.getCoarsening() == 0, "DEFAULT is not set.")

        self.assertTrue(so.getDiagnostics("num_iter") == None, "initial num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("num_inner_iter") == None, "initial num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("time") == None, "initial time is wrong.")
        self.assertTrue(so.getDiagnostics("set_up_time") == None, "initial set_up_time is wrong.")
        self.assertTrue(so.getDiagnostics("residual_norm") == None, "initial residual_norm is wrong.")
        self.assertTrue(so.getDiagnostics("converged") == None, "initial converged is wrong.")
        self.assertTrue(so.hasConverged() == None, "initial convergence flag is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_inner_iter") == 0, "initial cum_num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_iter") == 0, "initial cum_num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_time") ==0, "initial cum_time is wrong.")
        self.assertTrue(so.getDiagnostics("cum_set_up_time") == 0, "initial cum_set_up_time is wrong.")

        so._updateDiagnostics("num_iter",1)
        so._updateDiagnostics("num_inner_iter",2)
        so._updateDiagnostics("time",3)
        so._updateDiagnostics("set_up_time",4)
        so._updateDiagnostics("residual_norm",5)
        so._updateDiagnostics("converged",True)

        self.assertTrue(so.getDiagnostics("num_iter") == 1, "num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("num_inner_iter") == 2, "num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("time") == 3, "time is wrong.")
        self.assertTrue(so.getDiagnostics("set_up_time") == 4, "set_up_time is wrong.")
        self.assertTrue(so.getDiagnostics("residual_norm") == 5, "residual_norm is wrong.")
        self.assertTrue(so.getDiagnostics("converged"), "converged is wrong.")
        self.assertTrue(so.hasConverged(), "convergence flag is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_inner_iter") == 2, "cum_num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_iter") == 1, "cum_num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_time") ==3, "cum_time is wrong.")
        self.assertTrue(so.getDiagnostics("cum_set_up_time") == 4, "cum_set_up_time is wrong.")  
        
        so.resetDiagnostics()
        self.assertTrue(so.getDiagnostics("num_iter") == None, "initial num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("num_inner_iter") == None, "initial num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("time") == None, "initial time is wrong.")
        self.assertTrue(so.getDiagnostics("set_up_time") == None, "initial set_up_time is wrong.")
        self.assertTrue(so.getDiagnostics("residual_norm") == None, "initial residual_norm is wrong.")
        self.assertTrue(so.getDiagnostics("converged") == None, "initial converged is wrong.")
        self.assertTrue(so.hasConverged() == None, "initial convergence flag is wrong")       
        self.assertTrue(so.getDiagnostics("cum_num_inner_iter") == 2, "cum_num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_iter") == 1, "cum_num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_time") ==3, "cum_time is wrong.")
        self.assertTrue(so.getDiagnostics("cum_set_up_time") == 4, "cum_set_up_time is wrong.")

        so._updateDiagnostics("num_iter",10)
        so._updateDiagnostics("num_inner_iter",20)
        so._updateDiagnostics("time",30)
        so._updateDiagnostics("set_up_time",40)
        so._updateDiagnostics("residual_norm",50)
        so._updateDiagnostics("converged",False)

        self.assertTrue(so.getDiagnostics("num_iter") == 10, "num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("num_inner_iter") == 20, "num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("time") == 30, "time is wrong.")
        self.assertTrue(so.getDiagnostics("set_up_time") == 40, "set_up_time is wrong.")
        self.assertTrue(so.getDiagnostics("residual_norm") == 50, "residual_norm is wrong.")
        self.assertTrue(not so.getDiagnostics("converged"), "converged is wrong.")
        self.assertTrue(not so.hasConverged(), "convergence flag is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_inner_iter") == 22, "cum_num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_iter") == 11, "cum_num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_time") ==33, "cum_time is wrong.")
        self.assertTrue(so.getDiagnostics("cum_set_up_time") == 44, "cum_set_up_time is wrong.")  

        so.resetDiagnostics(all=True)
        self.assertTrue(so.getDiagnostics("num_iter") == None, "initial num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("num_inner_iter") == None, "initial num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("time") == None, "initial time is wrong.")
        self.assertTrue(so.getDiagnostics("set_up_time") == None, "initial set_up_time is wrong.")
        self.assertTrue(so.getDiagnostics("residual_norm") == None, "initial residual_norm is wrong.")
        self.assertTrue(so.getDiagnostics("converged") == None, "initial converged is wrong.")
        self.assertTrue(so.hasConverged() == None, "initial convergence flag is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_inner_iter") == 0, "initial cum_num_inner_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_num_iter") == 0, "initial cum_num_iter is wrong.")
        self.assertTrue(so.getDiagnostics("cum_time") ==0, "initial cum_time is wrong.")
        self.assertTrue(so.getDiagnostics("cum_set_up_time") == 0, "initial cum_set_up_time is wrong.")
        
    def test_setCoefficient_WithIllegalFunctionSpace(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        self.assertRaises(IllegalCoefficientFunctionSpace, mypde.setValue, C=Vector(0.,FunctionOnBoundary(self.domain)))

    def test_setCoefficient_WithWrongName(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        self.assertRaises(IllegalCoefficient, mypde.setValue, ROMA=0.)

    def test_resetCoefficient_WithWrongShape(self):
        mypde=LinearPDE(self.domain,numEquations=2,debug=self.DEBUG)
        self.assertRaises(IllegalCoefficientValue, mypde.setValue, C=0.)

    def test_reducedOn(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setReducedOrderOn()
        mypde.setValue(A=kronecker(self.domain),D=x[0],Y=x[0])
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')

    def test_attemptToChangeOrderAfterDefinedCoefficient(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        self.assertRaises(RuntimeError, mypde.setReducedOrderOn)

    def test_reducedOnConfig(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        self.assertEqual((mypde.getFunctionSpaceForSolution(), mypde.getFunctionSpaceForEquation()),(ReducedSolution(self.domain),ReducedSolution(self.domain)),"reduced function spaces expected.")
    #
    #  set coefficients for scalars:
    #
    def test_setCoefficient_A_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),Function(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A").isEmpty(),"A is empty after reset of right hand side coefficients")

    def test_setCoefficient_B_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B").isEmpty(),"B is empty after reset of right hand side coefficients")

    def test_setCoefficient_C_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C").isEmpty(),"C is empty after reset of right hand side coefficients")

    def test_setCoefficient_D_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        coeff=mypde.getCoefficient("D")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),Function(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D").isEmpty(),"D is empty after reset of right hand side coefficients")

    def test_setCoefficient_X_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),Function(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty after reset of right hand side coefficients")

    def test_setCoefficient_Y_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=1.)
        coeff=mypde.getCoefficient("Y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),Function(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty after reset of right hand side coefficients")

    def test_setCoefficient_y_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=1.)
        coeff=mypde.getCoefficient("y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty after reset of right hand side coefficients")

    def test_setCoefficient_d_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=1.)
        coeff=mypde.getCoefficient("d")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d").isEmpty(),"d is empty after reset of right hand side coefficients")

    def test_setCoefficient_d_contact_Scalar(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(d_contact=1.)
            coeff=mypde.getCoefficient("d_contact")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1,1))

            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is empty after reset of right hand side coefficients")

    def test_setCoefficient_y_contact_Scalar(self):
        d=self.domain.getDim()
        if self.domain.supportsContactElements():
            mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
            mypde.setValue(y_contact=1.)
            coeff=mypde.getCoefficient("y_contact")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1))

            mypde.resetRightHandSideCoefficients()
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty after reset of right hand side coefficients")

    def test_setCoefficient_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is empty after reset of right hand side coefficients")

    def test_setCoefficient_B_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is empty after reset of right hand side coefficients")

    def test_setCoefficient_C_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is empty after reset of right hand side coefficients")

    def test_setCoefficient_D_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=1.)
        coeff=mypde.getCoefficient("D_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is empty after reset of right hand side coefficients")

    def test_setCoefficient_X_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty after reset of right hand side coefficients")

    def test_setCoefficient_Y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=1.)
        coeff=mypde.getCoefficient("Y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty after reset of right hand side coefficients")

    def test_setCoefficient_y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=1.)
        coeff=mypde.getCoefficient("y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty after reset of right hand side coefficients")

    def test_setCoefficient_d_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=1.)
        coeff=mypde.getCoefficient("d_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty after reset of right hand side coefficients")

    def test_setCoefficient_d_contact_reduced_Scalar(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(d_contact_reduced=1.)
            coeff=mypde.getCoefficient("d_contact_reduced")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))

            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_reduced_Scalar(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
            mypde.setValue(y_contact_reduced=1.)
            coeff=mypde.getCoefficient("y_contact_reduced")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))

            mypde.resetRightHandSideCoefficients()
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty after reset of right hand side coefficients")

    def test_setCoefficient_r_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty after reset of right hand side coefficients")

    def test_setCoefficient_q_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")

    def test_setCoefficient_r_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty after reset of right hand side coefficients")

    def test_setCoefficient_q_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))

        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")

    def test_setCoefficient_A_reduced_Scalar_usingA(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((d,d)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='A'
            FS=Function
        else:
            coeff_name='A_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),FS(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_B_reduced_Scalar_usingB(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='B'
            FS=Function
        else:
            coeff_name='B_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),FS(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_C_reduced_Scalar_usingC(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='C'
            FS=Function
        else:
            coeff_name='C_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),FS(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_D_reduced_Scalar_usingD(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Scalar(1.,ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='D'
            FS=Function
        else:
            coeff_name='D_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FS(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_X_reduced_Scalar_usingX(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient('X_reduced')
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient('X_reduced').isEmpty(),"X_reduced is not empty after reset of right hand side coefficients")

    def test_setCoefficient_Y_reduced_Scalar_usingY(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Scalar(1.,ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='Y'
            FS=Function
        else:
            coeff_name='Y_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FS(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient(coeff_name).isEmpty(),"%s is not empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_y_reduced_Scalar_using_y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='y'
            FS=FunctionOnBoundary
        else:
            coeff_name='y_reduced'
            FS=ReducedFunctionOnBoundary
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FS(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient(coeff_name).isEmpty(),"%s is not empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_d_reduced_Scalar_using_d(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='d'
            FS=FunctionOnBoundary
        else:
            coeff_name='d_reduced'
            FS=ReducedFunctionOnBoundary
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FS(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_d_contact_reduced_Scalar_using_d_contact(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(d_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
            if self.specialInterpolationSupported():
                coeff_name='d_contact'
                FS=FunctionOnContactZero
            else:
                coeff_name='d_contact_reduced'
                FS=ReducedFunctionOnContactZero
            coeff=mypde.getCoefficient(coeff_name)
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FS(self.domain),1,1))
            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_y_contact_reduced_Scalar_using_y_contact(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
            mypde.setValue(y_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
            if self.specialInterpolationSupported():
                coeff_name='y_contact'
                FS=FunctionOnContactZero
            else:
                coeff_name='y_contact_reduced'
                FS=ReducedFunctionOnContactZero
            coeff=mypde.getCoefficient(coeff_name)
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FS(self.domain),1))
            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    #
    #  set coefficients for systems:
    #
    def test_setCoefficient_A_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A").isEmpty(),"A is empty after reset of right hand side coefficients")
        
    def test_setCoefficient_B_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B").isEmpty(),"B is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C").isEmpty(),"C is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D").isEmpty(),"D is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),Function(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),Function(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnBoundary(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty after reset of right hand side coefficients")
    def test_setCoefficient_d_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d").isEmpty(),"d is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_System(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(d_contact=numpy.ones((self.N,self.N)))
            coeff=mypde.getCoefficient("d_contact")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnContactZero(self.domain),self.N,self.N))
            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_System(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
            mypde.setValue(y_contact=numpy.ones((self.N,)))
            coeff=mypde.getCoefficient("y_contact")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnContactZero(self.domain),self.N))
            mypde.resetRightHandSideCoefficients()
            self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty after reset of right hand side coefficients")
    def test_setCoefficient_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_B_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_d_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_reduced_System(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(d_contact_reduced=numpy.ones((self.N,self.N)))
            coeff=mypde.getCoefficient("d_contact_reduced")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_reduced_System(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
            mypde.setValue(y_contact_reduced=numpy.ones((self.N,)))
            coeff=mypde.getCoefficient("y_contact_reduced")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
            mypde.resetRightHandSideCoefficients()
            self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_r_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty after reset of right hand side coefficients")
    def test_setCoefficient_q_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")
    def test_setCoefficient_r_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty after reset of right hand side coefficients")
    def test_setCoefficient_q_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")

    def test_setCoefficient_A_reduced_System_using_A(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((self.N,d,self.N,d)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='A'
            FS=Function
        else:
            coeff_name='A_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),FS(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_B_reduced_System_using_B(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((self.N,d,self.N)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='B'
            FS=Function
        else:
            coeff_name='B_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),FS(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_C_reduced_System_using_C(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((self.N,self.N,d)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='C'
            FS=Function
        else:
            coeff_name='C_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),FS(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_D_reduced_System_using_D(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Data(numpy.ones((self.N,self.N)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='D'
            FS=Function
        else:
            coeff_name='D_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FS(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_X_reduced_System_using_X(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=Data(numpy.ones((self.N,d)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='X'
            FS=Function
        else:
            coeff_name='X_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),FS(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient(coeff_name).isEmpty(),"%s is not empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_Y_reduced_System_using_Y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Data(numpy.ones((self.N,)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='Y'
            FS=Function
        else:
            coeff_name='Y_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FS(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient(coeff_name).isEmpty(),"%s is not empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_y_reduced_System_using_y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Data(numpy.ones((self.N,)),ReducedFunctionOnBoundary(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='y'
            FS=FunctionOnBoundary
        else:
            coeff_name='y_reduced'
            FS=ReducedFunctionOnBoundary
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FS(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient(coeff_name).isEmpty(),"%s is not empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_d_reduced_System_using_d(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='d'
            FS=FunctionOnBoundary
        else:
            coeff_name='d_reduced'
            FS=ReducedFunctionOnBoundary
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),FS(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_d_contact_reduced_System_using_d_contact(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(d_contact=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnContactZero(self.domain)))
            if self.specialInterpolationSupported():
                coeff_name='d_contact'
                FS=FunctionOnContactZero
            else:
                coeff_name='d_contact_reduced'
                FS=ReducedFunctionOnContactZero
            coeff=mypde.getCoefficient(coeff_name)
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),FS(self.domain),self.N,self.N))
            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_y_contact_reduced_System_using_y_contact(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
            mypde.setValue(y_contact=Data(numpy.ones((self.N,)),ReducedFunctionOnContactZero(self.domain)))
            if self.specialInterpolationSupported():
                coeff_name='y_contact'
                FS=FunctionOnContactZero
            else:
                coeff_name='y_contact_reduced'
                FS=ReducedFunctionOnContactZero
            coeff=mypde.getCoefficient(coeff_name)
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FS(self.domain),self.N))
            mypde.resetRightHandSideCoefficients()
            self.assertTrue(mypde.getCoefficient(coeff_name).isEmpty(),"%s is not empty after reset of right hand side coefficients"%coeff_name)

    def test_resetCoefficient_HomogeneousConstraint(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(A=kronecker(self.domain),Y=1.,q=whereZero(x[0]))
        u1=mypde.getSolution()
        mypde.setValue(Y=2.)
        u2=mypde.getSolution()
        self.assertTrue(self.check(u2,2*u1),'solution is wrong.')

    def test_resetCoefficient_InHomogeneousConstraint(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setSymmetryOn()
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        x=self.domain.getX()
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.,r=1,q=whereZero(x[0]))
        u1=mypde.getSolution()
        mypde.setValue(Y=2.,D=2)
        u2=mypde.getSolution()
        self.assertTrue(self.check(u2,u1),'first solution is wrong.')
        u2=mypde.getSolution()
        self.assertTrue(self.check(u2,u1),'first solution is wrong.')
        mypde.setValue(r=2,Y=4.)
        u2=mypde.getSolution()
        self.assertTrue(self.check(u2,2*u1),'second solution is wrong.')

    def test_Status(self):
        DIM=self.domain.getDim()
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSymmetryOn()
        mypde.getSolverOptions().setTolerance(self.RES_TOL)
        mypde.setValue(A=kronecker(self.domain), q=whereZero(x[0])+whereZero(x[0]-1.), Y=2.)
        x1=self.domain.getX()
        u1_ref=x1[0]*(1.-x1[0])
        u1=mypde.getSolution()
        error1=Lsup(u1-u1_ref)/Lsup(u1_ref)
        self.assertTrue(mypde.getDomainStatus() == mypde.getSystemStatus(), "status of first pde does not match domain status.")
        try:
            self.domain.setX(x*5)
        except:
            # setX not supported
            return
        self.assertTrue(mypde.getDomainStatus() != mypde.getSystemStatus(), "status of first pde matches updated domain status.")
        x2=self.domain.getX()
        u2_ref=x2[0]*(5.-x2[0])
        u2=mypde.getSolution()
        error2=Lsup(u2-u2_ref)/Lsup(u2_ref)
        self.assertTrue(error2 <= max(error1,self.RES_TOL)*10., "solution of second PDE wrong.")
        self.assertTrue(mypde.getDomainStatus() == mypde.getSystemStatus(), "status of second pde does not match domain status.")

    def test_symmetryCheckTrue_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        D=3*numpy.ones((self.N,self.N))
        d=4*numpy.ones((self.N,self.N))
        d_contact=5*numpy.ones((self.N,self.N))
        pars={"A":A, "B":B, "C":C, "D":D, "d":d, "A_reduced":-A, "B_reduced":-B, "C_reduced":-C, "D_reduced":-D, "d_reduced":-d}
        if self.domain.supportsContactElements():
                pars["d_contact"]=d_contact
                pars["d_contact_reduced"]=-d_contact
        mypde.setValue(**pars)
        self.assertTrue(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B=B,C=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D=D)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d=d)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_System(self):
        if self.domain.supportsContactElements():
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            d_contact=5*numpy.ones((self.N,self.N))
            d_contact[0,1]=0.
            mypde.setValue(d_contact=d_contact)
            self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A_reduced=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D_reduced=D)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d_reduced=d)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_reduced_System(self):
        if self.domain.supportsContactElements():
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            d_contact=5*numpy.ones((self.N,self.N))
            d_contact[0,1]=0.
            mypde.setValue(d_contact_reduced=d_contact)
            self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckTrue_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        D=3
        d=4
        d_contact=5
        pars={"A":A, "B":B, "C":C, "D":D, "d":d, "A_reduced":-A, "B_reduced":-B, "C_reduced":-C, "D_reduced":-D, "d_reduced":-d}
        if self.domain.supportsContactElements():
                pars["d_contact"]=d_contact
                pars["d_contact_reduced"]=-d_contact
        mypde.setValue(**pars)
        self.assertTrue(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B=B,C=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A_reduced=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    #
    #   solver checks (single PDE)
    #
    def test_symmetryOnIterative(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_symmetryOnDirect(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_DIRECT(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 	  
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_MINRES_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_MINRES_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_MINRES_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG',0):
                print("AMG test disabled on MPI build")
                return                
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_MINRES_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_MINRES_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_MINRES_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(50)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(50)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 	  
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            mypde.getSolverOptions().setTruncation(50)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(50)                         
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(50)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(50)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 	  
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')        
    def test_GMRES_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_GAUSS_SEIDEL(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_AMG(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 	  
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            mypde.getSolverOptions().setTruncation(10)
            mypde.getSolverOptions().setRestart(20)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    #
    #   solver checks (PDE system)
    #
    def test_symmetryOnIterative_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_symmetryOnDirect_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_JACOBI_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_GAUSS_SEIDEL_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_AMG_System(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return
            A=Tensor4(0.,Function(self.domain))
            D=Tensor(1.,Function(self.domain))
            Y=Vector(self.domain.getDim(),Function(self.domain))
            for i in range(self.domain.getDim()): 
                A[i,:,i,:]=kronecker(self.domain)
                D[i,i]+=i
                Y[i]+=i
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=A,D=D,Y=Y)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PCG_ILU0_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_DIRECT_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_JACOBI_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_GAUSS_SEIDEL_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_AMG_System(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 	  
            A=Tensor4(0.,Function(self.domain))
            D=Tensor(1.,Function(self.domain))
            Y=Vector(self.domain.getDim(),Function(self.domain))
            for i in range(self.domain.getDim()): 
                A[i,:,i,:]=kronecker(self.domain)
                D[i,i]+=i
                Y[i]+=i
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=A,D=D,Y=Y)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_ILU0_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_JACOBI_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_GAUSS_SEIDEL_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_AMG_System(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return
            A=Tensor4(0.,Function(self.domain))
            D=Tensor(1.,Function(self.domain))
            Y=Vector(self.domain.getDim(),Function(self.domain))
            for i in range(self.domain.getDim()): 
                A[i,:,i,:]=kronecker(self.domain)
                D[i,i]+=i
                Y[i]+=i
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=A,D=D,Y=Y)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_PRES20_ILU0_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_JACOBI_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        # u=mypde.getSolution(verbose=self.VERBOSE,truncation=5)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_GAUSS_SEIDEL_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        # u=mypde.getSolution(verbose=self.VERBOSE,truncation=5)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_AMG_System(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG',0):
                print("AMG test disabled on MPI build")
                return   	  
            A=Tensor4(0.,Function(self.domain))
            D=Tensor(1.,Function(self.domain))
            Y=Vector(self.domain.getDim(),Function(self.domain))
            for i in range(self.domain.getDim()): 
                A[i,:,i,:]=kronecker(self.domain)
                D[i,i]+=i
                Y[i]+=i
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=A,D=D,Y=Y)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            # u=mypde.getSolution(verbose=self.VERBOSE,truncation=5)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_ILU0_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        # u=mypde.getSolution(verbose=self.VERBOSE,truncation=5)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_JACOBI_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_GAUSS_SEIDEL_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_AMG_System(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 	  
            A=Tensor4(0.,Function(self.domain))
            D=Tensor(1.,Function(self.domain))
            Y=Vector(self.domain.getDim(),Function(self.domain))
            for i in range(self.domain.getDim()): 
                A[i,:,i,:]=kronecker(self.domain)
                D[i,i]+=i
                Y[i]+=i
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=A,D=D,Y=Y)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_ILU0_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_JACOBI_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_GAUSS_SEIDEL_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.GAUSS_SEIDEL)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_AMG_System(self):
        if self.order!=2:
            if getEscriptParamInt('DISABLE_AMG', 0):
                print("AMG test disabled on MPI build")
                return 	  
            A=Tensor4(0.,Function(self.domain))
            D=Tensor(1.,Function(self.domain))
            Y=Vector(self.domain.getDim(),Function(self.domain))
            for i in range(self.domain.getDim()): 
                A[i,:,i,:]=kronecker(self.domain)
                D[i,i]+=i
                Y[i]+=i
            mypde=LinearPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(A=A,D=D,Y=Y)
            mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
            mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
            mypde.getSolverOptions().setVerbosity(self.VERBOSE)
            mypde.getSolverOptions().setTruncation(10)
            mypde.getSolverOptions().setRestart(20)
            u=mypde.getSolution()
            self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_ILU0_System(self):
        A=Tensor4(0.,Function(self.domain))
        D=Tensor(1.,Function(self.domain))
        Y=Vector(self.domain.getDim(),Function(self.domain))
        for i in range(self.domain.getDim()): 
            A[i,:,i,:]=kronecker(self.domain)
            D[i,i]+=i
            Y[i]+=i
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=A,D=D,Y=Y)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
        mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        mypde.getSolverOptions().setTruncation(10)
        mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')


    def test_FluxScalar0(self):
        pde= LinearPDE(self.domain, numEquations=1, numSolutions=1)
        u=self.domain.getX()[0]
        f = pde.getFlux(u)
        self.assertEqual(f.getShape(),(self.domain.getDim(),),"wrong shape of result.")
        self.assertEqual(f.getFunctionSpace(),Function(self.domain),"wrong function space")
        self.assertEqual(Lsup(f),0.,"wrong result")

    def test_FluxScalar(self):
        pde= LinearPDE(self.domain, numEquations=1, numSolutions=1)
        pde.setValue(X=kronecker(self.domain)[0]*1., B=kronecker(self.domain)[1]*2, A=5*kronecker(self.domain))
        x=self.domain.getX()[0]
        f = pde.getFlux(x)
        self.assertEqual(f.getShape(),(self.domain.getDim(),),"wrong shape of result.")
        self.assertEqual(f.getFunctionSpace(),Function(self.domain),"wrong function space")
        f_ref=x*kronecker(self.domain)[1]*2+(5-1)*kronecker(self.domain)[0]
        self.assertTrue(self.check(f, f_ref),"wrong result")

    def test_FluxScalarReduced(self):
        pde= LinearPDE(self.domain, numEquations=1, numSolutions=1)
        pde.setValue(X_reduced=kronecker(self.domain)[0]*1., B_reduced=kronecker(self.domain)[1]*2, A_reduced=5*kronecker(self.domain))
        x=self.domain.getX()[0]
        f = pde.getFlux(x)
        self.assertEqual(f.getShape(),(self.domain.getDim(),),"wrong shape of result.")
        if self.specialInterpolationSupported():
            FS=Function
        else:
            FS=ReducedFunction
        self.assertEqual(f.getFunctionSpace(),FS(self.domain),"wrong function space")
        f_ref=Data(x*kronecker(self.domain)[1]*2+(5-1)*kronecker(self.domain)[0], ReducedFunction(self.domain))
        self.assertTrue(self.check(f, f_ref),"wrong result")

    def test_FluxSystem0(self):
        pde= LinearPDE(self.domain, numEquations=2, numSolutions=2)
        u=self.domain.getX()
        f = pde.getFlux(u)
        self.assertEqual(f.getShape(),(2, self.domain.getDim()),"wrong shape of result.")
        self.assertEqual(f.getFunctionSpace(),Function(self.domain),"wrong function space")
        self.assertEqual(Lsup(f),0.,"wrong result")

    def test_FluxSystem(self):
        pde= LinearPDE(self.domain, numEquations=2, numSolutions=2)
        X=Data(0., (2, self.domain.getDim()), Function(self.domain))
        X[0,0]=1
        B=Data(0., (2, self.domain.getDim(),2), Function(self.domain))
        B[0,0,0]=5
        A=Data(0., (2, self.domain.getDim(),2, self.domain.getDim()), Function(self.domain))
        A[0,0,0,0]=10
        pde.setValue(X=X, B=B, A=A)
        x=self.domain.getX()
        f = pde.getFlux(x[:2])
        self.assertEqual(f.getShape(),(2, self.domain.getDim()),"wrong shape of result.")
        self.assertEqual(f.getFunctionSpace(),Function(self.domain),"wrong function space")
        f_ref=X*(5*x[0]-1+10)
        self.assertTrue(self.check(f, f_ref),"wrong result")
    def test_FluxSystemReduced(self):
        pde= LinearPDE(self.domain, numEquations=2, numSolutions=2)
        X=Data(0., (2, self.domain.getDim()), ReducedFunction(self.domain))
        X[0,0]=1
        B=Data(0., (2, self.domain.getDim(),2), ReducedFunction(self.domain))
        B[0,0,0]=5
        A=Data(0., (2, self.domain.getDim(),2, self.domain.getDim()), ReducedFunction(self.domain))
        A[0,0,0,0]=10
        pde.setValue(X=X, B=B, A=A)
        x=self.domain.getX()
        f = pde.getFlux(x[:2])
        self.assertEqual(f.getShape(),(2, self.domain.getDim()),"wrong shape of result.")
        if self.specialInterpolationSupported():
            FS=Function
        else:
            FS=ReducedFunction
        self.assertEqual(f.getFunctionSpace(),FS(self.domain),"wrong function space")
        f_ref=X*(5*x[0]-1+10)
        self.assertTrue(self.check(f, f_ref),"wrong result")
        
class Test_LinearPDE(Test_LinearPDE_noLumping):
    def test_Lumping_attemptToSetA(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
           mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
           mypde.setValue(A=kronecker(self.domain))
           mypde.getSolverOptions().setVerbosity(self.VERBOSE)
           u=mypde.getSolution()        
        except ValueError:
           success=False
        self.assertTrue(not success,'error should be issued')
    def test_Lumping_attemptToSetB(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
           mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
           mypde.setValue(B=kronecker(self.domain)[0])
           mypde.getSolverOptions().setVerbosity(self.VERBOSE)
           u=mypde.getSolution()
        except ValueError:
           success=False
        self.assertTrue(not success,'error should be issued')
    def test_Lumping_attemptToSetC(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
           mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
           mypde.setValue(C=kronecker(self.domain)[0])
           mypde.getSolverOptions().setVerbosity(self.VERBOSE)
           u=mypde.getSolution()
        except ValueError:
           success=False
        self.assertTrue(not success,'error should be issued')
        
    def test_Lumping_attemptToSetA_reduced(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
           mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
           mypde.setValue(A_reduced=kronecker(self.domain))
           mypde.getSolverOptions().setVerbosity(self.VERBOSE)
           u=mypde.getSolution()
        except ValueError:
           success=False
        self.assertTrue(not success,'error should be issued')
    def test_Lumping_attemptToSetB_reduced(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
           mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
           mypde.setValue(B_reduced=kronecker(self.domain)[0])
           mypde.getSolverOptions().setVerbosity(self.VERBOSE)
           u=mypde.getSolution()
        except ValueError:
           success=False
        self.assertTrue(not success,'error should be issued')
    def test_Lumping_attemptToSetC_reduced(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
           mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
           mypde.setValue(C_reduced=kronecker(self.domain)[0])
           mypde.getSolverOptions().setVerbosity(self.VERBOSE)
           u=mypde.getSolution()
        except ValueError:
           success=False
        self.assertTrue(not success,'error should be issued')
        
    def test_Lumping(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')
    def test_Constrained_Lumping(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=1.,Y=1.,q=whereZero(x[0]),r=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'solution is wrong.')

    def test_Lumping_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=numpy.array([[1.,0.],[0.,2.]]),Y=numpy.array([1.,2.]))
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,numpy.ones((2,))),'solution is wrong.')
    def test_Constrained_Lumping_System(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=numpy.array([[1.,0.],[0.,2.]]),Y=numpy.array([1.,2.]), \
                       q=whereZero(x[0])*[0.,1],r=[0.,1.])
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,numpy.ones((2,))),'solution is wrong.')

    def test_Lumping_updateRHS(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,1.),'first solution is wrong.')
        mypde.setValue(Y=2.,q=whereZero(x[0]),r=2.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,2.),'second solution is wrong.')
    def test_Lumping_updateOperator(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        mypde.setValue(D=2.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.assertTrue(self.check(u,0.5),'second solution is wrong.')


class Test_TransportPDE(Test_linearPDEs):
    N=4

    def test_init_Init(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)

    def test_setCoefficient_WithWrongName(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        self.assertRaises(IllegalCoefficient, mypde.setValue, ROMA=Vector(0.,FunctionOnBoundary(self.domain)))

    def test_setCoefficient_WithIllegalFunctionSpace(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        self.assertRaises(IllegalCoefficientFunctionSpace, mypde.setValue,C=Vector(0.,FunctionOnBoundary(self.domain)))
        
    def test_resetCoefficient_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numEquations=2,debug=self.DEBUG)
        self.assertRaises(IllegalCoefficientValue, mypde.setValue, C=0.)

    def test_setInitialSolution_scalar(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        mypde.setInitialSolution(1.)

    def test_setInitialSolution_scalar_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        self.assertRaises(ValueError,mypde.setInitialSolution,[1.,2.])

    def test_setInitialSolution_system(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        mypde.setInitialSolution([1.,2.])

    def test_setInitialSolution_system_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        self.assertRaises(ValueError, mypde.setInitialSolution,1.)

    def test_attemptToChangeOrderAfterDefinedCoefficient(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        self.assertRaises(RuntimeError, mypde.setReducedOrderOn)

    def test_reducedOnConfig(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        self.assertEqual((mypde.getFunctionSpaceForSolution(), mypde.getFunctionSpaceForEquation()),(ReducedSolution(self.domain),ReducedSolution(self.domain)),"reduced function spaces expected.")
    #
    #  set coefficients for scalars:
    #
    def test_setCoefficient_M_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=1.)
        coeff=mypde.getCoefficient("M")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),Function(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("M").isEmpty(),"M is empty after reset of right hand side coefficients")

    def test_setCoefficient_A_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),Function(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A").isEmpty(),"A is empty after reset of right hand side coefficients")
    def test_setCoefficient_B_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B").isEmpty(),"B is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C").isEmpty(),"C is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        coeff=mypde.getCoefficient("D")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),Function(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D").isEmpty(),"D is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),Function(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=1.)
        coeff=mypde.getCoefficient("Y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),Function(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=1.)
        coeff=mypde.getCoefficient("y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty after reset of right hand side coefficients")
    def test_setCoefficient_d_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=1.)
        coeff=mypde.getCoefficient("d")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d").isEmpty(),"d is empty after reset of right hand side coefficients")
    def test_setCoefficient_m_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=1.)
        coeff=mypde.getCoefficient("m")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("m").isEmpty(),"m is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_Scalar(self):
        if self.domain.supportsContactElements():
            d=self.domain.getDim()
            mypde=TransportPDE(self.domain,debug=self.DEBUG)
            mypde.setValue(d_contact=1.)
            coeff=mypde.getCoefficient("d_contact")
            self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1,1))
            mypde.resetRightHandSideCoefficients()
            self.assertFalse(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=1.)
        coeff=mypde.getCoefficient("y_contact")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty after reset of right hand side coefficients")

    def test_setCoefficient_M_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M_reduced=1.)
        coeff=mypde.getCoefficient("M_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("M_reduced").isEmpty(),"M_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_B_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=1.)
        coeff=mypde.getCoefficient("D_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=1.)
        coeff=mypde.getCoefficient("Y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=1.)
        coeff=mypde.getCoefficient("y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_m_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m_reduced=1.)
        coeff=mypde.getCoefficient("m_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("m_reduced").isEmpty(),"m_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=1.)
        coeff=mypde.getCoefficient("d_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=1.)
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=1.)
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_r_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty after reset of right hand side coefficients")
    def test_setCoefficient_q_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")
    def test_setCoefficient_r_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is empty after reset of right hand side coefficients")
    def test_setCoefficient_q_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")

    def test_setCoefficient_M_reduced_Scalar_usingM(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=Scalar(1.,ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='M'
            FS=Function
        else:
            coeff_name='M_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FS(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_A_reduced_Scalar_usingA(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((d,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_B_reduced_Scalar_usingB(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_reduced_Scalar_usingC(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_reduced_Scalar_usingD(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_reduced_Scalar_usingX(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_reduced_Scalar_usingY(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_reduced_Scalar_using_y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_m_reduced_Scalar_using_m(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='m'
            FS=FunctionOnBoundary
        else:
            coeff_name='m_reduced'
            FS=ReducedFunctionOnBoundary
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FS(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_d_reduced_Scalar_using_d(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_reduced_Scalar_using_d_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_contact_reduced").isEmpty(),"M is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_reduced_Scalar_using_y_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty after reset of right hand side coefficients")
    #
    #  set coefficients for systems:
    #
    def test_setCoefficient_M_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("M")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("M").isEmpty(),"M is empty after reset of right hand side coefficients")
    def test_setCoefficient_A_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A").isEmpty(),"A is empty after reset of right hand side coefficients")
    def test_setCoefficient_B_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B").isEmpty(),"B is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C").isEmpty(),"C is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D").isEmpty(),"D is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),Function(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is not empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),Function(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y").isEmpty(),"Y is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnBoundary(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y").isEmpty(),"y is not empty after reset of right hand side coefficients")
    def test_setCoefficient_m_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("m")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("m").isEmpty(),"m is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d").isEmpty(),"d is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_System(self):
        
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnContactZero(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnContactZero(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty after reset of right hand side coefficients")
    def test_setCoefficient_M_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("M_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("M_reduced").isEmpty(),"M_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_B_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X").isEmpty(),"X is empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_m_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("m_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("m_reduced").isEmpty(),"m_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"X is not empty after reset of right hand side coefficients")
    def test_setCoefficient_r_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is not empty after reset of right hand side coefficients")
    def test_setCoefficient_q_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")
    def test_setCoefficient_r_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("r").isEmpty(),"r is no empty after reset of right hand side coefficients")
    def test_setCoefficient_q_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("q").isEmpty(),"q is empty after reset of right hand side coefficients")

    def test_setCoefficient_M_reduced_System_using_M(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=Data(numpy.ones((self.N,self.N)),ReducedFunction(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='M'
            FS=Function
        else:
            coeff_name='M_reduced'
            FS=ReducedFunction
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FS(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_A_reduced_System_using_A(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((self.N,d,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_B_reduced_System_using_B(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((self.N,d,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_C_reduced_System_using_C(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((self.N,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_D_reduced_System_using_D(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Data(numpy.ones((self.N,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_X_reduced_System_using_X(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=Data(numpy.ones((self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_Y_reduced_System_using_Y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Data(numpy.ones((self.N,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_y_reduced_System_using_y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Data(numpy.ones((self.N,)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty after reset of right hand side coefficients")
    def test_setCoefficient_m_reduced_System_using_m(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        if self.specialInterpolationSupported():
            coeff_name='m'
            FS=FunctionOnBoundary
        else:
            coeff_name='m_reduced'
            FS=ReducedFunctionOnBoundary
        coeff=mypde.getCoefficient(coeff_name)
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FS(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient(coeff_name).isEmpty(),"%s is empty after reset of right hand side coefficients"%coeff_name)

    def test_setCoefficient_d_reduced_System_using_d(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_d_contact_reduced_System_using_d_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertFalse(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is empty after reset of right hand side coefficients")
    def test_setCoefficient_y_contact_reduced_System_using_y_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Data(numpy.ones((self.N,)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.assertEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
        mypde.resetRightHandSideCoefficients()
        self.assertTrue(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty after reset of right hand side coefficients")

    def test_symmetryCheckTrue_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=100*numpy.ones((self.N,self.N))
        A=numpy.ones((self.N,d,self.N,d))
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        D=3*numpy.ones((self.N,self.N))
        d=4*numpy.ones((self.N,self.N))
        m=64*numpy.ones((self.N,self.N))
        d_contact=5*numpy.ones((self.N,self.N))
        mypde.setValue(M=M,A=A,B=B,C=C,D=D,d=d,d_contact=d_contact,m=m,M_reduced=-M,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact, m_reduced=-m)
        self.assertTrue(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_M_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=numpy.ones((self.N,self.N))
        M[1,0]=0.
        mypde.setValue(M=M)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_BC_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B=B,C=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D=D)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_m_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        m=4*numpy.ones((self.N,self.N))
        m[0,1]=0.
        mypde.setValue(m=m)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d=d)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numpy.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact=d_contact)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_M_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=3*numpy.ones((self.N,self.N))
        M[0,1]=0.
        mypde.setValue(M_reduced=M)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A_reduced=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_BC_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D_reduced=D)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_m_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        m=4*numpy.ones((self.N,self.N))
        m[0,1]=0.
        mypde.setValue(m_reduced=m)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d_reduced=d)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numpy.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact_reduced=d_contact)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    #==============================================================
    def test_symmetryCheckTrue_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=100
        A=numpy.ones((d,d))
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        D=3
        m=10
        d=4
        d_contact=5
        mypde.setValue(M=M,A=A,B=B,C=C,D=D,d=d,m=m,d_contact=d_contact,M_reduced=-M,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact,m_reduced=-m)
        self.assertTrue(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B=B,C=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A_reduced=A)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.assertTrue(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_reducedOn(self):
        dt=0.1
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setInitialSolution(10.)
        mypde.setValue(M=1.,Y=1)
        u=mypde.getSolution(dt)
        self.assertTrue(u.getFunctionSpace() == ReducedSolution(self.domain), "wrong function space")
        self.assertTrue(self.check(u,10.+dt),'solution is wrong.')

    def Off_test_reducedOff(self):
        dt=0.1
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        mypde.setInitialSolution(10.)
        mypde.setValue(M=1.,Y=1.)
        u=mypde.getSolution(0.1)
        self.assertTrue(u.getFunctionSpace() == Solution(self.domain), "wrong function space")
        self.assertTrue(self.check(u,10.+dt),'solution is wrong.')
