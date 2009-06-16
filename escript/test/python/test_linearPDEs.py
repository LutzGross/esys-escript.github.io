
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
Test suite for linearPDEs class

The tests must be linked with a Domain class object in the setUp method:

   from esys.finley import Rectangle
   class Test_LinearPDEOnFinley(Test_LinearPDE):
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

from esys.escript.util import Lsup,kronecker,interpolate,whereZero, outer, swap_axes
from esys.escript import Function,FunctionOnBoundary,FunctionOnContactZero,Solution,ReducedSolution,Vector,ContinuousFunction,Scalar, ReducedFunction,ReducedFunctionOnBoundary,ReducedFunctionOnContactZero,Data, Tensor4, Tensor
from esys.escript.linearPDEs import LinearPDE,IllegalCoefficientValue,Poisson, IllegalCoefficientFunctionSpace, TransportPDE, IllegalCoefficient, Helmholtz, LameEquation, SolverOptions
import numpy
import unittest

class Test_linearPDEs(unittest.TestCase):
    TOL=1.e-6
    SOLVER_TOL=1.e-10
    DEBUG=False
    VERBOSE=False
    def check(self,arg,ref_arg,tol=None):
        """
        checks if arg and ref_arg are nearly identical using the L{Lsup<esys.escript.util.Lsup>}
        """
        if tol==None: tol=self.TOL
        return Lsup(arg-ref_arg)<=tol*Lsup(ref_arg)
    
class Test_LameEquation(Test_linearPDEs):

    def test_config(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        d=self.domain.getDim()
        self.failUnlessEqual((mypde.getNumEquations(), mypde.getNumSolutions(), mypde.getSolverOptions().isSymmetric()),(d,d,True),"set up incorrect")

    def test_setCoefficient_q(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(q=x)

        q_ref=interpolate(x,Solution(self.domain))
        self.failUnless(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(self.check(mypde.getCoefficient("q"),q_ref),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_r(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(r=x)

        r_ref=interpolate(x,Solution(self.domain))
        self.failUnless(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(self.check(mypde.getCoefficient("r"),r_ref),"r is nor x")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")


    def test_setCoefficient_F(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(F=x)

        Y_ref=interpolate(x,Function(self.domain))
        self.failUnless(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(self.check(mypde.getCoefficient("Y"),Y_ref),"Y is not x")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_f(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(f=x)

        y_ref=interpolate(x,FunctionOnBoundary(self.domain))
        self.failUnless(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(self.check(mypde.getCoefficient("y"),y_ref),"d is not x[0]")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_sigma(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(sigma=outer(x,x))

        X_ref=interpolate(outer(x,x),Function(self.domain))
        self.failUnless(self.check(mypde.getCoefficient("A"),0),"A is not 0")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(self.check(mypde.getCoefficient("X"),X_ref),"X is not x X x")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_lambda(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(lame_lambda=x[0])


        k3=kronecker(Function(self.domain))
        k3Xk3=outer(k3,k3)
        A_ref=x[0]*k3Xk3

        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_mu(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(lame_mu=x[0])


        k3=kronecker(Function(self.domain))
        k3Xk3=outer(k3,k3)
        A_ref=x[0]*(swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3))

        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_lambdamu(self):
        mypde=LameEquation(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(lame_lambda=x[0], lame_mu=x[1])

        k3=kronecker(Function(self.domain))
        k3Xk3=outer(k3,k3)
        A_ref=x[0]*k3Xk3+x[1]*(swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3))

        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

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
       self.failUnless(self.check(u,u_ex,10*self.TOL),"incorrect solution")

class Test_Helmholtz(Test_linearPDEs):

    def test_config(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        self.failUnlessEqual((mypde.getNumEquations(), mypde.getNumSolutions(), mypde.getSolverOptions().isSymmetric()),(1,1,True),"set up incorrect")
    def test_setCoefficient_q(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(q=whereZero(x[0]))

        q_ref=interpolate(whereZero(x[0]),Solution(self.domain))
        A_ref=kronecker(self.domain)

        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(self.check(mypde.getCoefficient("q"),q_ref),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_r(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(r=x[0])

        r_ref=interpolate(x[0],Solution(self.domain))
        A_ref=kronecker(self.domain)
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(self.check(mypde.getCoefficient("r"),r_ref),"r is nor x[0]")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")


    def test_setCoefficient_f(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(f=x[0])

        Y_ref=interpolate(x[0],Function(self.domain))
        A_ref=kronecker(self.domain)
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(self.check(mypde.getCoefficient("Y"),Y_ref),"Y is not x[0]")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_alpha(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(alpha=x[0])

        d_ref=interpolate(x[0],FunctionOnBoundary(self.domain))
        A_ref=kronecker(self.domain)
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(self.check(mypde.getCoefficient("d"),d_ref),"d is not x[0]")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_g(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(g=x[0])

        y_ref=interpolate(x[0],FunctionOnBoundary(self.domain))
        A_ref=kronecker(self.domain)
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(self.check(mypde.getCoefficient("y"),y_ref),"y is not x[0]")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_setCoefficient_omega(self):
        mypde=Helmholtz(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(omega=x[0])

        D_ref=interpolate(x[0],Function(self.domain))
        A_ref=kronecker(self.domain)
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(self.check(mypde.getCoefficient("D"),D_ref),"D is not x[0]")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")

    def test_solve(self):
       d=self.domain.getDim()
       cf=ContinuousFunction(self.domain)
       u_ex=Scalar(1.,cf)
       mypde=Helmholtz(self.domain)
       mypde.setValue(f=3,omega=3,alpha=2,g=2)
       u=mypde.getSolution()
       self.failUnless(self.check(u,u_ex,10*self.TOL),"incorrect solution")

class Test_Poisson(Test_linearPDEs):

    def test_config(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        self.failUnlessEqual((mypde.getNumEquations(), mypde.getNumSolutions(), mypde.getSolverOptions().isSymmetric()),(1,1,True),"set up incorrect")
    def test_setCoefficient_q(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        q_ref=interpolate(whereZero(x[0]),Solution(self.domain))
        A_ref=kronecker(self.domain)
        mypde.setValue(q=whereZero(x[0]))
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(self.check(mypde.getCoefficient("q"),q_ref),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")
    def test_setCoefficient_f(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        Y_ref=interpolate(x[0],Function(self.domain))
        A_ref=kronecker(self.domain)
        mypde.setValue(f=x[0])
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(self.check(mypde.getCoefficient("Y"),Y_ref),"Y is not x[0]")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(mypde.getCoefficient("Y_reduced").isEmpty(),"Y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")
    def test_setCoefficient_f_reduced(self):
        mypde=Poisson(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        Y_ref=interpolate(x[0],ReducedFunction(self.domain))
        A_ref=kronecker(self.domain)
        mypde.setValue(f_reduced=x[0])
        self.failUnless(self.check(mypde.getCoefficient("A"),A_ref),"A is not kronecker")
        self.failUnless(mypde.getCoefficient("B").isEmpty(),"B is not empty")
        self.failUnless(mypde.getCoefficient("C").isEmpty(),"C is not empty")
        self.failUnless(mypde.getCoefficient("D").isEmpty(),"D is not empty")
        self.failUnless(mypde.getCoefficient("X").isEmpty(),"X is not empty")
        self.failUnless(mypde.getCoefficient("Y").isEmpty(),"Y is not empty")
        self.failUnless(mypde.getCoefficient("y").isEmpty(),"y is not empty")
        self.failUnless(mypde.getCoefficient("d").isEmpty(),"d is not empty")
        self.failUnless(mypde.getCoefficient("d_contact").isEmpty(),"d_contact is not empty")
        self.failUnless(mypde.getCoefficient("y_contact").isEmpty(),"y_contact is not empty")
        self.failUnless(mypde.getCoefficient("A_reduced").isEmpty(),"A_reduced is not empty")
        self.failUnless(mypde.getCoefficient("B_reduced").isEmpty(),"B_reduced is not empty")
        self.failUnless(mypde.getCoefficient("C_reduced").isEmpty(),"C_reduced is not empty")
        self.failUnless(mypde.getCoefficient("D_reduced").isEmpty(),"D_reduced is not empty")
        self.failUnless(mypde.getCoefficient("X_reduced").isEmpty(),"X_reduced is not empty")
        self.failUnless(self.check(mypde.getCoefficient("Y_reduced"),Y_ref),"Y_reduced is not x[0]")
        self.failUnless(mypde.getCoefficient("y_reduced").isEmpty(),"y_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_reduced").isEmpty(),"d_reduced is not empty")
        self.failUnless(mypde.getCoefficient("d_contact_reduced").isEmpty(),"d_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("y_contact_reduced").isEmpty(),"y_contact_reduced is not empty")
        self.failUnless(mypde.getCoefficient("q").isEmpty(),"q is not empty")
        self.failUnless(mypde.getCoefficient("r").isEmpty(),"r is not empty")
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
       self.failUnless(self.check(u,u_ex,10*self.TOL),"incorrect solution")

class Test_LinearPDE_noLumping(Test_linearPDEs):
    N=4
    def test_SolverOptions(self):
        so=SolverOptions()

        self.failUnless(so.getLevelMax() == 10, "initial  LevelMax is wrong.")
        self.failUnlessRaises(ValueError,so.setLevelMax,-1)
        so.setLevelMax(3)
        self.failUnless(so.getLevelMax() == 3, "LevelMax is wrong.")

        self.failUnless(so.getCoarseningThreshold() == 0.05, "initial  CoarseningThreshold is wrong.")
        self.failUnlessRaises(ValueError,so.setCoarseningThreshold,-1)
        so.setCoarseningThreshold(0.1)
        self.failUnless(so.getCoarseningThreshold() == 0.1, "CoarseningThreshold is wrong.")

        self.failUnless(so.getNumSweeps() == 2, "initial  Sweeps is wrong.")
        self.failUnlessRaises(ValueError,so.setNumSweeps,-1)
        so.setNumSweeps(3)
        self.failUnless(so.getNumSweeps() == 3, "Sweeps is wrong.")

        self.failUnless(so.getNumPreSweeps() == 2, "initial  PreSweeps is wrong.")
        self.failUnlessRaises(ValueError,so.setNumPreSweeps,-1)
        so.setNumPreSweeps(4)
        self.failUnless(so.getNumPreSweeps() == 4, "PreSweeps is wrong.")

        self.failUnless(so.getNumPostSweeps() == 2, "initial  PreSweeps is wrong.")
        self.failUnlessRaises(ValueError,so.setNumPostSweeps,-1)
        so.setNumPostSweeps(5)
        self.failUnless(so.getNumPostSweeps() == 5, "PreSweeps is wrong.")

        self.failUnless(so.getTolerance() == 1.e-8, "initial Tolerance is wrong.")
        self.failUnlessRaises(ValueError,so.setTolerance,-1)
        so.setTolerance(0.2)
        self.failUnless(so.getTolerance() == 0.2, "Tolerance is wrong.")

        self.failUnless(so.getAbsoluteTolerance() == 0., "initial  AbsoluteTolerance is wrong.")
        self.failUnlessRaises(ValueError,so.setAbsoluteTolerance,-1)
        so.setAbsoluteTolerance(0.3)
        self.failUnless(so.getAbsoluteTolerance() == 0.3, "AbsoluteTolerance is wrong.")

        self.failUnless(so.getInnerTolerance() == 0.9, "initial InnerTolerance is wrong.")
        self.failUnlessRaises(ValueError,so.setInnerTolerance,-1)
        so.setInnerTolerance(0.4)
        self.failUnless(so.getInnerTolerance() == 0.4, "InnerTolerance is wrong.")

        self.failUnless(so.getDropTolerance() == 0.01, "initial DropTolerance is wrong.")
        self.failUnlessRaises(ValueError,so.setDropTolerance,-1)
        so.setDropTolerance(0.5)
        self.failUnless(so.getDropTolerance() == 0.5, "DropDropTolerance is wrong.")

        self.failUnless(so.getDropStorage() == 2., "initial DropStorage is wrong.")
        self.failUnlessRaises(ValueError,so.setDropStorage,-1)
        so.setDropStorage(10)
        self.failUnless(so.getDropStorage() == 10, "DropStorage is wrong.")
        
        self.failUnless(so.getRelaxationFactor() == 0.3, "initial RelaxationFactor is wrong.")
        self.failUnlessRaises(ValueError,so.setRelaxationFactor,-1)
        so.setRelaxationFactor(0.1)
        self.failUnless(so.getRelaxationFactor() == 0.1, "Relaxation is wrong.")


        self.failUnless(so.getIterMax() == 10000, "initial IterMax is wrong.")
        self.failUnlessRaises(ValueError,so.setIterMax,0)
        so.setIterMax(11)
        self.failUnless(so.getIterMax() == 11, "IterMax is wrong.")

        self.failUnless(so.getInnerIterMax() == 10, "initial InnerIterMax is wrong.")
        self.failUnlessRaises(ValueError,so.setInnerIterMax,0)
        so.setInnerIterMax(12)
        self.failUnless(so.getInnerIterMax() == 12, "InnerIterMax is wrong.")

        self.failUnless(so.getTruncation() == 20, "initial Truncation is wrong.")
        self.failUnlessRaises(ValueError,so.setTruncation,0)
        so.setTruncation(13)
        self.failUnless(so.getTruncation() == 13, "Truncation is wrong.")

        self.failUnless(so.getRestart() == None, "initial Truncation is wrong.")
        self.failUnlessRaises(ValueError,so.setTruncation,0)
        so.setRestart(14)
        self.failUnless(so.getRestart() == 14, "Truncation is wrong.")
        so.setRestart(None)
        self.failUnless(so.getRestart() == None, "Truncation is wrong.")
	
        self.failUnless(not so.isVerbose(), "initial verbosity flag is wrong.")
        so.setVerbosityOn()
        self.failUnless(so.isVerbose(), "verbosity (1) flag is wrong.")
        so.setVerbosityOff()
        self.failUnless(not so.isVerbose(), "verbosity (2) flag is wrong.")
        so.setVerbosity(verbose=True)
        self.failUnless(so.isVerbose(), "verbosity (3) flag is wrong.")
        so.setVerbosity(verbose=False)
        self.failUnless(not so.isVerbose(), "verbosity (4) flag is wrong.")

        self.failUnless(not so.isSymmetric(), "initial symmetry flag is wrong.")
        so.setSymmetryOn()
        self.failUnless(so.isSymmetric(), "symmetry (1) flag is wrong.")
        so.setSymmetryOff()
        self.failUnless(not so.isSymmetric(), "symmetry (2) flag is wrong.")
        so.setSymmetry(flag=True)
        self.failUnless(so.isSymmetric(), "symmetry (3) flag is wrong.")
        so.setSymmetry(flag=False)
        self.failUnless(not so.isSymmetric(), "symmetry (4) flag is wrong.")

        self.failUnless(so.adaptInnerTolerance(), "initial InnerToleranceAdaption flag is wrong.")
        so.setInnerToleranceAdaptionOn()
        self.failUnless(so.adaptInnerTolerance(), "InnerToleranceAdaption (1) flag is wrong.")
        so.setInnerToleranceAdaptionOff()
        self.failUnless(not so.adaptInnerTolerance(), "InnerToleranceAdaption (2) flag is wrong.")
        so.setInnerToleranceAdaption(adapt=True)
        self.failUnless(so.adaptInnerTolerance(), "InnerToleranceAdaption (3) flag is wrong.")
        so.setInnerToleranceAdaption(adapt=False)
        self.failUnless(not so.adaptInnerTolerance(), "InnerToleranceAdaption (4) flag is wrong.")
     
        self.failUnless(not so.acceptConvergenceFailure(), "initial acceptConvergenceFailure flag is wrong.")
        so.setAcceptanceConvergenceFailureOn()
        self.failUnless(so.acceptConvergenceFailure(), "acceptConvergenceFailure (1) flag is wrong.")
        so.setAcceptanceConvergenceFailureOff()
        self.failUnless(not so.acceptConvergenceFailure(), "acceptConvergenceFailure (2) flag is wrong.")
        so.setAcceptanceConvergenceFailure(accept=True)
        self.failUnless(so.acceptConvergenceFailure(), "acceptConvergenceFailure (3) flag is wrong.")
        so.setAcceptanceConvergenceFailure(accept=False)
        self.failUnless(not so.acceptConvergenceFailure(), "acceptConvergenceFailure (4) flag is wrong.")   
        
        self.failUnless(so.getReordering() == 30, "initial Reordering is wrong.")
        self.failUnlessRaises(ValueError,so.setReordering,-1)
        so.setReordering(so.NO_REORDERING)
        self.failUnless(so.getReordering() == 17, "NO_REORDERING is not set.")
        so.setReordering(so.MINIMUM_FILL_IN)
        self.failUnless(so.getReordering() == 18, "MINIMUM_FILL_IN is not set.")
        so.setReordering(so.NESTED_DISSECTION)
        self.failUnless(so.getReordering() == 19, "NESTED_DISSECTION is not set.")
        so.setReordering(so.DEFAULT_REORDERING)
        self.failUnless(so.getReordering() == 30, "DEFAULT_REORDERING is not set.")
        
        self.failUnless(so.getPackage() == 0, "initial solver package is wrong.")
        self.failUnlessRaises(ValueError,so.setPackage,-1)
        so.setPackage(so.PASO)
        self.failUnless(so.getPackage() == 21, "PASO is not set.")
        so.setPackage(so.SUPER_LU)
        self.failUnless(so.getPackage() == 31, "SUPER_LU is not set.")
        so.setPackage(so.PASTIX)
        self.failUnless(so.getPackage() == 32, "PASTIX is not set.")
        so.setPackage(so.MKL)
        self.failUnless(so.getPackage() == 15, "MKL is not set.")
        so.setPackage(so.UMFPACK)
        self.failUnless(so.getPackage() == 16, "UMFPACK is not set.")
        so.setPackage(so.TRILINOS)
        self.failUnless(so.getPackage() == 24, "TRILINOS is not set.")

        self.failUnless(so.getSolverMethod() == 0, "initial SolverMethod is wrong.")
        self.failUnlessRaises(ValueError,so.setSolverMethod,-1)
        so.setSolverMethod(so.DIRECT)
        self.failUnless(so.getSolverMethod() == 1, "DIRECT is not set.")
        so.setSolverMethod(so.CHOLEVSKY)
        self.failUnless(so.getSolverMethod() == 2, "CHOLEVSKY is not set.")
        so.setSolverMethod(so.PCG)
        self.failUnless(so.getSolverMethod() == 3, "PCG is not set.")
        so.setSolverMethod(so.CR)
        self.failUnless(so.getSolverMethod() == 4, "CR is not set.")
        so.setSolverMethod(so.CGS)
        self.failUnless(so.getSolverMethod() == 5, "CGS is not set.")
        so.setSolverMethod(so.BICGSTAB)
        self.failUnless(so.getSolverMethod() == 6, "BICGSTAB is not set.")
        so.setSolverMethod(so.SSOR)
        self.failUnless(so.getSolverMethod() == 7, "SSOR is not set.")
        so.setSolverMethod(so.GMRES)
        self.failUnless(so.getSolverMethod() == 11, "GMRES is not set.")
        so.setSolverMethod(so.PRES20)
        self.failUnless(so.getSolverMethod() == 12, "PRES20 is not set.")
        so.setSolverMethod(so.LUMPING)
        self.failUnless(so.getSolverMethod() == 13, "LUMPING is not set.")
        so.setSolverMethod(so.ITERATIVE)
        self.failUnless(so.getSolverMethod() == 20, "ITERATIVE is not set.")
        so.setSolverMethod(so.AMG)
        self.failUnless(so.getSolverMethod() == 22, "AMG is not set.")
        so.setSolverMethod(so.NONLINEAR_GMRES)
        self.failUnless(so.getSolverMethod() == 25, "NONLINEAR_GMRES is not set.")
        so.setSolverMethod(so.TFQMR)
        self.failUnless(so.getSolverMethod() == 26, "TFQMR is not set.")
        so.setSolverMethod(so.MINRES)
        self.failUnless(so.getSolverMethod() == 27, "MINRES is not set.")
        so.setSolverMethod(so.GAUSS_SEIDEL)
        self.failUnless(so.getSolverMethod() == 28, "GAUSS_SEIDEL is not set.")
        so.setSolverMethod(so.DEFAULT)
        self.failUnless(so.getSolverMethod() == 0, "DEFAULT is not set.")

        self.failUnless(so.getPreconditioner() == 10, "initial Preconditioner is wrong.")
        self.failUnlessRaises(ValueError,so.setPreconditioner,-1)
        so.setPreconditioner(so.ILU0)
        self.failUnless(so.getPreconditioner() == 8, "ILU0 is not set.")
        so.setPreconditioner(so.SSOR)
        self.failUnless(so.getPreconditioner() == 7, "SSOR is not set.")
        so.setPreconditioner(so.ILUT)
        self.failUnless(so.getPreconditioner() == 9, "ILUT is not set.")
        so.setPreconditioner(so.JACOBI)
        self.failUnless(so.getPreconditioner() == 10, "JACOBI is not set.")
        so.setPreconditioner(so.AMG)
        self.failUnless(so.getPreconditioner() == 22, "AMG is not set.")
        so.setPreconditioner(so.REC_ILU)
        self.failUnless(so.getPreconditioner() == 23, "REC_ILU is not set.")
        so.setPreconditioner(so.GAUSS_SEIDEL)
        self.failUnless(so.getPreconditioner() == 28, "GAUSS_SEIDEL is not set.")
        so.setPreconditioner(so.RILU)
        self.failUnless(so.getPreconditioner() == 29, "RILU is not set.")
        so.setPreconditioner(so.NO_PRECONDITIONER)
        self.failUnless(so.getPreconditioner() == 36, "NO_PRECONDITIONER is not set.")        

        self.failUnless(so.getCoarsening() == 0, "initial Coarseningr is wrong.")
        self.failUnlessRaises(ValueError,so.setCoarsening,-1)
        so.setCoarsening(so.YAIR_SHAPIRA_COARSENING)
        self.failUnless(so.getCoarsening() == 33, "YAIR_SHAPIRA_COARSENING is not set.")
        so.setCoarsening(so.RUGE_STUEBEN_COARSENING)
        self.failUnless(so.getCoarsening() == 34, "RUGE_STUEBEN_COARSENING is not set.")
        so.setCoarsening(so.AGGREGATION_COARSENING)
        self.failUnless(so.getCoarsening() == 35, "AGREGATION_COARSENING is not set.")
        so.setCoarsening(so.DEFAULT)
        self.failUnless(so.getCoarsening() == 0, "DEFAULT is not set.")

        self.failUnless(so.getDiagnostics("num_iter") == None, "initial num_iter is wrong.")
        self.failUnless(so.getDiagnostics("num_inner_iter") == None, "initial num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("time") == None, "initial time is wrong.")
        self.failUnless(so.getDiagnostics("set_up_time") == None, "initial set_up_time is wrong.")
        self.failUnless(so.getDiagnostics("residual_norm") == None, "initial residual_norm is wrong.")
        self.failUnless(so.getDiagnostics("converged") == None, "initial converged is wrong.")
        self.failUnless(so.hasConverged() == None, "initial convergence flag is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_inner_iter") == 0, "initial cum_num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_iter") == 0, "initial cum_num_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_time") ==0, "initial cum_time is wrong.")
        self.failUnless(so.getDiagnostics("cum_set_up_time") == 0, "initial cum_set_up_time is wrong.")

        so._updateDiagnostics("num_iter",1)
        so._updateDiagnostics("num_inner_iter",2)
        so._updateDiagnostics("time",3)
        so._updateDiagnostics("set_up_time",4)
        so._updateDiagnostics("residual_norm",5)
        so._updateDiagnostics("converged",True)

        self.failUnless(so.getDiagnostics("num_iter") == 1, "num_iter is wrong.")
        self.failUnless(so.getDiagnostics("num_inner_iter") == 2, "num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("time") == 3, "time is wrong.")
        self.failUnless(so.getDiagnostics("set_up_time") == 4, "set_up_time is wrong.")
        self.failUnless(so.getDiagnostics("residual_norm") == 5, "residual_norm is wrong.")
        self.failUnless(so.getDiagnostics("converged"), "converged is wrong.")
        self.failUnless(so.hasConverged(), "convergence flag is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_inner_iter") == 2, "cum_num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_iter") == 1, "cum_num_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_time") ==3, "cum_time is wrong.")
        self.failUnless(so.getDiagnostics("cum_set_up_time") == 4, "cum_set_up_time is wrong.")  
        
        so.resetDiagnostics()
        self.failUnless(so.getDiagnostics("num_iter") == None, "initial num_iter is wrong.")
        self.failUnless(so.getDiagnostics("num_inner_iter") == None, "initial num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("time") == None, "initial time is wrong.")
        self.failUnless(so.getDiagnostics("set_up_time") == None, "initial set_up_time is wrong.")
        self.failUnless(so.getDiagnostics("residual_norm") == None, "initial residual_norm is wrong.")
        self.failUnless(so.getDiagnostics("converged") == None, "initial converged is wrong.")
        self.failUnless(so.hasConverged() == None, "initial convergence flag is wrong")       
        self.failUnless(so.getDiagnostics("cum_num_inner_iter") == 2, "cum_num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_iter") == 1, "cum_num_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_time") ==3, "cum_time is wrong.")
        self.failUnless(so.getDiagnostics("cum_set_up_time") == 4, "cum_set_up_time is wrong.")

        so._updateDiagnostics("num_iter",10)
        so._updateDiagnostics("num_inner_iter",20)
        so._updateDiagnostics("time",30)
        so._updateDiagnostics("set_up_time",40)
        so._updateDiagnostics("residual_norm",50)
        so._updateDiagnostics("converged",False)

        self.failUnless(so.getDiagnostics("num_iter") == 10, "num_iter is wrong.")
        self.failUnless(so.getDiagnostics("num_inner_iter") == 20, "num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("time") == 30, "time is wrong.")
        self.failUnless(so.getDiagnostics("set_up_time") == 40, "set_up_time is wrong.")
        self.failUnless(so.getDiagnostics("residual_norm") == 50, "residual_norm is wrong.")
        self.failUnless(not so.getDiagnostics("converged"), "converged is wrong.")
        self.failUnless(not so.hasConverged(), "convergence flag is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_inner_iter") == 22, "cum_num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_iter") == 11, "cum_num_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_time") ==33, "cum_time is wrong.")
        self.failUnless(so.getDiagnostics("cum_set_up_time") == 44, "cum_set_up_time is wrong.")  

        so.resetDiagnostics(all=True)
        self.failUnless(so.getDiagnostics("num_iter") == None, "initial num_iter is wrong.")
        self.failUnless(so.getDiagnostics("num_inner_iter") == None, "initial num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("time") == None, "initial time is wrong.")
        self.failUnless(so.getDiagnostics("set_up_time") == None, "initial set_up_time is wrong.")
        self.failUnless(so.getDiagnostics("residual_norm") == None, "initial residual_norm is wrong.")
        self.failUnless(so.getDiagnostics("converged") == None, "initial converged is wrong.")
        self.failUnless(so.hasConverged() == None, "initial convergence flag is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_inner_iter") == 0, "initial cum_num_inner_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_num_iter") == 0, "initial cum_num_iter is wrong.")
        self.failUnless(so.getDiagnostics("cum_time") ==0, "initial cum_time is wrong.")
        self.failUnless(so.getDiagnostics("cum_set_up_time") == 0, "initial cum_set_up_time is wrong.")
        
    def test_setCoefficient_WithIllegalFunctionSpace(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficientFunctionSpace, mypde.setValue, C=Vector(0.,FunctionOnBoundary(self.domain)))

    def test_setCoefficient_WithWrongName(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficient, mypde.setValue, ROMA=0.)

    def test_resetCoefficient_WithWrongShape(self):
        mypde=LinearPDE(self.domain,numEquations=2,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficientValue, mypde.setValue, C=0.)

    def test_reducedOn(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setReducedOrderOn()
        mypde.setValue(A=kronecker(self.domain),D=x[0],Y=x[0])
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')

    def test_attemptToChangeOrderAfterDefinedCoefficient(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        self.failUnlessRaises(RuntimeError, mypde.setReducedOrderOn)

    def test_reducedOnConfig(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        self.failUnlessEqual((mypde.getFunctionSpaceForSolution(), mypde.getFunctionSpaceForEquation()),(ReducedSolution(self.domain),ReducedSolution(self.domain)),"reduced function spaces expected.")
    #
    #  set coefficients for scalars:
    #
    def test_setCoefficient_A_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),Function(self.domain),1,1))
    def test_setCoefficient_B_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_C_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_D_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),Function(self.domain),1,1))
    def test_setCoefficient_X_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),Function(self.domain),1))
    def test_setCoefficient_Y_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=1.)
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),Function(self.domain),1))
    def test_setCoefficient_y_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=1.)
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=1.)
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=1.)
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=1.)
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1))
    def test_setCoefficient_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=1.)
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=1.)
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=1.)
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=1.)
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=1.)
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=1.)
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    def test_setCoefficient_r_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_q_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_r_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))
    def test_setCoefficient_q_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))

    def test_setCoefficient_A_reduced_Scalar_usingA(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((d,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar_usingB(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar_usingC(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar_usingD(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar_usingX(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar_usingY(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar_using_y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_reduced_Scalar_using_d(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar_using_d_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar_using_y_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    #
    #  set coefficients for systems:
    #
    def test_setCoefficient_A_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_B_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_C_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_D_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_X_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),Function(self.domain),self.N))
    def test_setCoefficient_Y_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),Function(self.domain),self.N))
    def test_setCoefficient_y_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_d_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_d_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_r_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_q_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_r_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
    def test_setCoefficient_q_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))

    def test_setCoefficient_A_reduced_System_using_A(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((self.N,d,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System_using_B(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((self.N,d,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System_using_C(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((self.N,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced_using_D(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Data(numpy.ones((self.N,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced_using_X(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=Data(numpy.ones((self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced_using_Y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Data(numpy.ones((self.N,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_reduced_System_using_y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Data(numpy.ones((self.N,)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_d_reduced_System_using_d(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System_using_d_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System_using_y_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Data(numpy.ones((self.N,)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
    def test_resetCoefficient_HomogeneousConstraint(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        x=self.domain.getX()
        mypde.setValue(A=kronecker(self.domain),Y=1.,q=whereZero(x[0]))
        u1=mypde.getSolution()
        mypde.setValue(Y=2.)
        u2=mypde.getSolution()
        self.failUnless(self.check(u2,2*u1),'solution is wrong.')

    def test_resetCoefficient_InHomogeneousConstraint(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setSymmetryOn()
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        x=self.domain.getX()
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.,r=1,q=whereZero(x[0]))
        u1=mypde.getSolution()
        mypde.setValue(Y=2.,D=2)
        u2=mypde.getSolution()
        self.failUnless(self.check(u2,u1),'first solution is wrong.')
        u2=mypde.getSolution()
        self.failUnless(self.check(u2,u1),'first solution is wrong.')
        mypde.setValue(r=2,Y=4.)
        u2=mypde.getSolution()
        self.failUnless(self.check(u2,2*u1),'second solution is wrong.')

    def test_symmetryCheckTrue_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        D=3*numpy.ones((self.N,self.N))
        d=4*numpy.ones((self.N,self.N))
        d_contact=5*numpy.ones((self.N,self.N))
        mypde.setValue(A=A,B=B,C=C,D=D,d=d,d_contact=d_contact,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact)
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numpy.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D_reduced=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d_reduced=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numpy.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact_reduced=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckTrue_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        D=3
        d=4
        d_contact=5
        mypde.setValue(A=A,B=B,C=C,D=D,d=d,d_contact=d_contact,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact)
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    #
    #   solver checks (single PDE)
    #
    def test_symmetryOnIterative(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_symmetryOnDirect(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PCG_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PCG_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PCG_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PCG_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_DIRECT(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_MINRES_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_MINRES_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_MINRES_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_MINRES_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PRES20_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PRES20_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PRES20_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PRES20_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.PRES20)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(50)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(50)                         
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(50)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(50)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(10)
	mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.ILU0)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(10)
	mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.RILU)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(10)
	mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_REC_ILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.GMRES)
	mypde.getSolverOptions().setPreconditioner(SolverOptions.REC_ILU)
	mypde.getSolverOptions().setVerbosity(self.VERBOSE)
	mypde.getSolverOptions().setTruncation(10)
	mypde.getSolverOptions().setRestart(20)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')
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
        self.failUnless(self.check(u,1.),'solution is wrong.')

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
        self.failUnless(not success,'error should be issued')
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
        self.failUnless(not success,'error should be issued')
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
        self.failUnless(not success,'error should be issued')
        
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
        self.failUnless(not success,'error should be issued')
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
        self.failUnless(not success,'error should be issued')
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
        self.failUnless(not success,'error should be issued')
        
    def test_Lumping(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_Constrained_Lumping(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=1.,Y=1.,q=whereZero(x[0]),r=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'solution is wrong.')

    def test_Lumping_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=numpy.array([[1.,0.],[0.,2.]]),Y=numpy.array([1.,2.]))
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,numpy.ones((2,))),'solution is wrong.')
    def test_Constrained_Lumping_System(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=numpy.array([[1.,0.],[0.,2.]]),Y=numpy.array([1.,2.]), \
                       q=whereZero(x[0])*[0.,1],r=[0.,1.])
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,numpy.ones((2,))),'solution is wrong.')

    def test_Lumping_updateRHS(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,1.),'first solution is wrong.')
        mypde.setValue(Y=2.,q=whereZero(x[0]),r=2.)
        mypde.getSolverOptions().setVerbosity(self.VERBOSE)
        u=mypde.getSolution()
        self.failUnless(self.check(u,2.),'second solution is wrong.')
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
        self.failUnless(self.check(u,0.5),'second solution is wrong.')


class Test_TransportPDE(Test_linearPDEs):
    N=4
    def test_init_useBackwardEuler(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG, useBackwardEuler=True)
        self.failUnless(mypde.useBackwardEuler()==True,'backward Euler should be used')
    def test_init_donntUseBackwardEuler(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG, useBackwardEuler=False)
        self.failUnless(mypde.useBackwardEuler()==False,'backward Euler should not be used')
    def test_setCoefficient_WithWrongName(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficient)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.setValue, ROMA=Vector(0.,FunctionOnBoundary(self.domain)))

    def test_setCoefficient_WithIllegalFunctionSpace(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficientFunctionSpace)
	mypde.getSolverOptions().setSolverMethod(SolverOptions.setValue,C=Vector(0.,FunctionOnBoundary(self.domain)))
        
    def test_resetCoefficient_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numEquations=2,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficientValue, mypde.setValue, C=0.)

    def test_setInitialSolution_scalar(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        mypde.setInitialSolution(1.)

    def test_setInitialSolution_scalar_negative(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        self.failUnlessRaises(RuntimeError, mypde.setInitialSolution,-1.)

    def test_setInitialSolution_scalar_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        self.failUnlessRaises(ValueError,mypde.setInitialSolution,[1.,2.])

    def test_setInitialSolution_system(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        mypde.setInitialSolution([1.,2.])

    def test_setInitialSolution_system(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        self.failUnlessRaises(RuntimeError, mypde.setInitialSolution,[-1,2.])

    def test_setInitialSolution_system_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        self.failUnlessRaises(ValueError, mypde.setInitialSolution,1.)


    def test_attemptToChangeOrderAfterDefinedCoefficient(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        self.failUnlessRaises(RuntimeError, mypde.setReducedOrderOn)

    def test_reducedOnConfig(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        self.failUnlessEqual((mypde.getFunctionSpaceForSolution(), mypde.getFunctionSpaceForEquation()),(ReducedSolution(self.domain),ReducedSolution(self.domain)),"reduced function spaces expected.")
    #
    #  set coefficients for scalars:
    #
    def test_setCoefficient_M_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=1.)
        coeff=mypde.getCoefficient("M")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),Function(self.domain),1,1))
    def test_setCoefficient_A_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),Function(self.domain),1,1))
    def test_setCoefficient_B_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_C_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_D_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),Function(self.domain),1,1))
    def test_setCoefficient_X_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),Function(self.domain),1))
    def test_setCoefficient_Y_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=1.)
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),Function(self.domain),1))
    def test_setCoefficient_y_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=1.)
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=1.)
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_m_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=1.)
        coeff=mypde.getCoefficient("m")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=1.)
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=1.)
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1))

    def test_setCoefficient_M_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M_reduced=1.)
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((d,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=1.)
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((d,)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=1.)
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=1.)
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_m_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m_reduced=1.)
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=1.)
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=1.)
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=1.)
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    def test_setCoefficient_r_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_q_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_r_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))
    def test_setCoefficient_q_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))

    def test_setCoefficient_M_reduced_Scalar_usingM(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_A_reduced_Scalar_usingA(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((d,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar_usingB(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar_usingC(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar_usingD(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar_usingX(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=Data(numpy.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar_usingY(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar_using_y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_m_reduced_Scalar_using_m(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_reduced_Scalar_using_d(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar_using_d_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar_using_y_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    #
    #  set coefficients for systems:
    #
    def test_setCoefficient_M_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("M")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_A_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_B_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_C_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_D_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_X_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),Function(self.domain),self.N))
    def test_setCoefficient_Y_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),Function(self.domain),self.N))
    def test_setCoefficient_y_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_m_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("m")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),FunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),FunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_M_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numpy.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numpy.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numpy.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numpy.ones((self.N,d)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_m_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=numpy.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_r_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_q_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_r_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
    def test_setCoefficient_q_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=numpy.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))

    def test_setCoefficient_M_System_reduced_using_D(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=Data(numpy.ones((self.N,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_A_reduced_System_using_A(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numpy.ones((self.N,d,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System_using_B(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numpy.ones((self.N,d,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System_using_C(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numpy.ones((self.N,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced_using_D(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Data(numpy.ones((self.N,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced_using_X(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=Data(numpy.ones((self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced_using_Y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Data(numpy.ones((self.N,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_reduced_System_using_y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Data(numpy.ones((self.N,)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_m_reduced_System_using_m(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_reduced_System_using_d(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System_using_d_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Data(numpy.ones((self.N,self.N)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumSolutions(), mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System_using_y_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Data(numpy.ones((self.N,)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(), mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))

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
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_M_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=numpy.ones((self.N,self.N))
        M[1,0]=0.
        mypde.setValue(M=M)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_BC_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_m_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        m=4*numpy.ones((self.N,self.N))
        m[0,1]=0.
        mypde.setValue(m=m)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numpy.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_M_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=3*numpy.ones((self.N,self.N))
        M[0,1]=0.
        mypde.setValue(M_reduced=M)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_BC_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((self.N,self.N,d))
        B=2*numpy.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        D=3*numpy.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D_reduced=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_m_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        m=4*numpy.ones((self.N,self.N))
        m[0,1]=0.
        mypde.setValue(m_reduced=m)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d=4*numpy.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d_reduced=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numpy.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact_reduced=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

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
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numpy.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numpy.ones((d,))
        B=2*numpy.ones((d,))
        B[0]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_reducedOn(self):
        dt=0.1
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setInitialSolution(10.)
        mypde.setValue(M=1.,Y=1)
        u=mypde.getSolution(dt)
        self.failUnless(u.getFunctionSpace() == ReducedSolution(self.domain), "wrong function space")
        self.failUnless(self.check(u,10.+dt),'solution is wrong.')

    def Off_test_reducedOff(self):
        dt=0.1
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        mypde.setInitialSolution(10.)
        mypde.setValue(M=1.,Y=1.)
        u=mypde.getSolution(0.1)
        self.failUnless(u.getFunctionSpace() == Solution(self.domain), "wrong function space")
        self.failUnless(self.check(u,10.+dt),'solution is wrong.')
