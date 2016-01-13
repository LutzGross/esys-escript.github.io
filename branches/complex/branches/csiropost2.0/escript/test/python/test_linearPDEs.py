
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
from esys.escript.linearPDEs import LinearPDE,IllegalCoefficientValue,Poisson, IllegalCoefficientFunctionSpace, TransportPDE, IllegalCoefficient, Helmholtz, LameEquation
import numarray
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
        self.failUnlessEqual((mypde.getNumEquations(),mypde.getNumSolutions(),mypde.isSymmetric()),(d,d,True),"set up incorrect")

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
        self.failUnlessEqual((mypde.getNumEquations(),mypde.getNumSolutions(),mypde.isSymmetric()),(1,1,True),"set up incorrect")
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
        self.failUnlessEqual((mypde.getNumEquations(),mypde.getNumSolutions(),mypde.isSymmetric()),(1,1,True),"set up incorrect")
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
    def test_setCoefficient_WithIllegalFunctionSpace(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficientFunctionSpace,mypde.setValue, C=Vector(0.,FunctionOnBoundary(self.domain)))

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
        self.failUnlessRaises(RuntimeError,mypde.setReducedOrderOn)

    def test_reducedOnConfig(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        self.failUnlessEqual((mypde.getFunctionSpaceForSolution(),mypde.getFunctionSpaceForEquation()),(ReducedSolution(self.domain),ReducedSolution(self.domain)),"reduced function spaces expected.")
    #
    #  set coefficients for scalars:
    #
    def test_setCoefficient_A_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numarray.ones((d,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,d),Function(self.domain),1,1))
    def test_setCoefficient_B_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numarray.ones((d,)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_C_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numarray.ones((d,)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_D_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),Function(self.domain),1,1))
    def test_setCoefficient_X_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numarray.ones((d,)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((d,),Function(self.domain),1))
    def test_setCoefficient_Y_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=1.)
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),Function(self.domain),1))
    def test_setCoefficient_y_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=1.)
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=1.)
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=1.)
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=1.)
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1))
    def test_setCoefficient_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numarray.ones((d,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numarray.ones((d,)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numarray.ones((d,)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=1.)
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numarray.ones((d,)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=1.)
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=1.)
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=1.)
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=1.)
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=1.)
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    def test_setCoefficient_r_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_q_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_r_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))
    def test_setCoefficient_q_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))

    def test_setCoefficient_A_reduced_Scalar_usingA(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numarray.ones((d,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar_usingB(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numarray.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar_usingC(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numarray.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar_usingD(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar_usingX(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=Data(numarray.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar_usingY(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar_using_y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_reduced_Scalar_using_d(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar_using_d_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar_using_y_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    #
    #  set coefficients for systems:
    #
    def test_setCoefficient_A_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numarray.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_B_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numarray.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_C_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numarray.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_D_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_X_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numarray.ones((self.N,d)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,d),Function(self.domain),self.N))
    def test_setCoefficient_Y_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),Function(self.domain),self.N))
    def test_setCoefficient_y_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),FunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_d_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),FunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),FunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numarray.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numarray.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numarray.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numarray.ones((self.N,d)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_System_reduced(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_d_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_r_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_q_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_r_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
    def test_setCoefficient_q_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))

    def test_setCoefficient_A_reduced_System_using_A(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numarray.ones((self.N,d,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System_using_B(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numarray.ones((self.N,d,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System_using_C(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numarray.ones((self.N,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced_using_D(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Data(numarray.ones((self.N,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced_using_X(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=Data(numarray.ones((self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced_using_Y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Data(numarray.ones((self.N,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_reduced_System_using_y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Data(numarray.ones((self.N,)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_d_reduced_System_using_d(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Data(numarray.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System_using_d_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Data(numarray.ones((self.N,self.N)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System_using_y_contact(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Data(numarray.ones((self.N,)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
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
        x=self.domain.getX()
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.,r=1,q=whereZero(x[0]))
        u1=mypde.getSolution(verbose=self.VERBOSE)
        mypde.setValue(Y=2.,D=2)
        u2=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u2,u1),'first solution is wrong.')
        u2=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u2,u1),'first solution is wrong.')
        mypde.setValue(r=2,Y=4.)
        u2=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u2,2*u1),'second solution is wrong.')

    def test_symmetryCheckTrue_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((self.N,d,self.N,d))
        C=2*numarray.ones((self.N,self.N,d))
        B=2*numarray.ones((self.N,d,self.N))
        D=3*numarray.ones((self.N,self.N))
        d=4*numarray.ones((self.N,self.N))
        d_contact=5*numarray.ones((self.N,self.N))
        mypde.setValue(A=A,B=B,C=C,D=D,d=d,d_contact=d_contact,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact)
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((self.N,self.N,d))
        B=2*numarray.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        D=3*numarray.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d=4*numarray.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numarray.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_System(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((self.N,self.N,d))
        B=2*numarray.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        D=3*numarray.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D_reduced=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d=4*numarray.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d_reduced=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_reduced_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numarray.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact_reduced=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckTrue_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((d,d))
        C=2*numarray.ones((d,))
        B=2*numarray.ones((d,))
        D=3
        d=4
        d_contact=5
        mypde.setValue(A=A,B=B,C=C,D=D,d=d,d_contact=d_contact,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact)
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((d,))
        B=2*numarray.ones((d,))
        B[0]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((d,))
        B=2*numarray.ones((d,))
        B[0]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    #
    #   solver checks (single PDE)
    #
    def test_symmetryOnIterative(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_symmetryOnDirect(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.setSolverMethod(mypde.DIRECT)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PCG_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.setSolverMethod(mypde.PCG,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PCG_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.setSolverMethod(mypde.PCG,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PCG_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.setSolverMethod(mypde.PCG,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_DIRECT(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        mypde.setSolverMethod(mypde.DIRECT)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.BICGSTAB,mypde.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.BICGSTAB,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_BICGSTAB_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.BICGSTAB,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_MINRES_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.MINRES,mypde.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_MINRES_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.MINRES,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_MINRES_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.MINRES,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.TFQMR,mypde.JACOBI)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.TFQMR,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_TFQMR_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.TFQMR,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PRES20_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.PRES20,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PRES20_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.PRES20,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_PRES20_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.PRES20,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=50)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=50)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRESnoRestart_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=50)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_JACOBI(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=10,restart=20)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_ILU0(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=10,restart=20)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_GMRES_truncation_restart_RILU(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=kronecker(self.domain),D=1.,Y=1.)
	mypde.setSolverMethod(mypde.GMRES,mypde.RILU)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=10,restart=20)
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
        u=mypde.getSolution(verbose=self.VERBOSE)
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
        mypde.setSolverMethod(mypde.DIRECT)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
        mypde.setSolverMethod(mypde.PCG,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
        mypde.setSolverMethod(mypde.PCG,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
        mypde.setSolverMethod(mypde.DIRECT)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.BICGSTAB,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.BICGSTAB,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.PRES20,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.PRES20,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.GMRES,mypde.JACOBI)
        # u=mypde.getSolution(verbose=self.VERBOSE,truncation=5)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.GMRES,mypde.ILU0)
        # u=mypde.getSolution(verbose=self.VERBOSE,truncation=5)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.GMRES,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.GMRES,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
	mypde.setSolverMethod(mypde.GMRES,mypde.JACOBI)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=10,restart=20)
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
	mypde.setSolverMethod(mypde.GMRES,mypde.ILU0)
        u=mypde.getSolution(verbose=self.VERBOSE,truncation=10,restart=20)
        self.failUnless(self.check(u,1.),'solution is wrong.')

class Test_LinearPDE(Test_LinearPDE_noLumping):
    def test_Lumping_attemptToSetA(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
	   mypde.setSolverMethod(mypde.LUMPING)
           mypde.setValue(A=kronecker(self.domain))
           u=mypde.getSolution(verbose=self.VERBOSE)
        except ValueError:
           success=False
        self.failUnless(not success,'error should be issued')
    def test_Lumping_attemptToSetB(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
	   mypde.setSolverMethod(mypde.LUMPING)
           mypde.setValue(B=kronecker(self.domain)[0])
           u=mypde.getSolution(verbose=self.VERBOSE)
        except ValueError:
           success=False
        self.failUnless(not success,'error should be issued')
    def test_Lumping_attemptToSetC(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
	   mypde.setSolverMethod(mypde.LUMPING)
           mypde.setValue(C=kronecker(self.domain)[0])
           u=mypde.getSolution(verbose=self.VERBOSE)
        except ValueError:
           success=False
        self.failUnless(not success,'error should be issued')
        
    def test_Lumping_attemptToSetA_reduced(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
	   mypde.setSolverMethod(mypde.LUMPING)
           mypde.setValue(A_reduced=kronecker(self.domain))
           u=mypde.getSolution(verbose=self.VERBOSE)
        except ValueError:
           success=False
        self.failUnless(not success,'error should be issued')
    def test_Lumping_attemptToSetB_reduced(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
	   mypde.setSolverMethod(mypde.LUMPING)
           mypde.setValue(B_reduced=kronecker(self.domain)[0])
           u=mypde.getSolution(verbose=self.VERBOSE)
        except ValueError:
           success=False
        self.failUnless(not success,'error should be issued')
    def test_Lumping_attemptToSetC_reduced(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
        try: 
           success=True
	   mypde.setSolverMethod(mypde.LUMPING)
           mypde.setValue(C_reduced=kronecker(self.domain)[0])
           u=mypde.getSolution(verbose=self.VERBOSE)
        except ValueError:
           success=False
        self.failUnless(not success,'error should be issued')
        
    def test_Lumping(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')
    def test_Constrained_Lumping(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.LUMPING)
        mypde.setValue(D=1.,Y=1.,q=whereZero(x[0]),r=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'solution is wrong.')

    def test_Lumping_System(self):
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.LUMPING)
        mypde.setValue(D=numarray.array([[1.,0.],[0.,2.]]),Y=numarray.array([1.,2.]))
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,numarray.ones((2,))),'solution is wrong.')
    def test_Constrained_Lumping_System(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.LUMPING)
        mypde.setValue(D=numarray.array([[1.,0.],[0.,2.]]),Y=numarray.array([1.,2.]), \
                       q=whereZero(x[0])*[0.,1],r=[0.,1.])
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,numarray.ones((2,))),'solution is wrong.')

    def test_Lumping_updateRHS(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,1.),'first solution is wrong.')
        mypde.setValue(Y=2.,q=whereZero(x[0]),r=2.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        self.failUnless(self.check(u,2.),'second solution is wrong.')
    def test_Lumping_updateOperator(self):
        x=self.domain.getX()
        mypde=LinearPDE(self.domain,debug=self.DEBUG)
	mypde.setSolverMethod(mypde.LUMPING)
        mypde.setValue(D=1.,Y=1.)
        u=mypde.getSolution(verbose=self.VERBOSE)
        mypde.setValue(D=2.)
        u=mypde.getSolution(verbose=self.VERBOSE)
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
        self.failUnlessRaises(IllegalCoefficient,mypde.setValue, ROMA=Vector(0.,FunctionOnBoundary(self.domain)))

    def test_setCoefficient_WithIllegalFunctionSpace(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficientFunctionSpace,mypde.setValue,C=Vector(0.,FunctionOnBoundary(self.domain)))
        
    def test_resetCoefficient_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numEquations=2,debug=self.DEBUG)
        self.failUnlessRaises(IllegalCoefficientValue, mypde.setValue, C=0.)

    def test_setInitialSolution_scalar(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        mypde.setInitialSolution(1.)

    def test_setInitialSolution_scalar_negative(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        self.failUnlessRaises(RuntimeError,mypde.setInitialSolution,-1.)

    def test_setInitialSolution_scalar_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numSolutions=1,debug=self.DEBUG)
        self.failUnlessRaises(ValueError,mypde.setInitialSolution,[1.,2.])

    def test_setInitialSolution_system(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        mypde.setInitialSolution([1.,2.])

    def test_setInitialSolution_system(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        self.failUnlessRaises(RuntimeError,mypde.setInitialSolution,[-1,2.])

    def test_setInitialSolution_system_WithWrongShape(self):
        mypde=TransportPDE(self.domain,numSolutions=2,debug=self.DEBUG)
        self.failUnlessRaises(ValueError,mypde.setInitialSolution,1.)


    def test_attemptToChangeOrderAfterDefinedCoefficient(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        self.failUnlessRaises(RuntimeError, mypde.setReducedOrderOn)

    def test_reducedOnConfig(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        self.failUnlessEqual((mypde.getFunctionSpaceForSolution(),mypde.getFunctionSpaceForEquation()),(ReducedSolution(self.domain),ReducedSolution(self.domain)),"reduced function spaces expected.")
    #
    #  set coefficients for scalars:
    #
    def test_setCoefficient_M_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=1.)
        coeff=mypde.getCoefficient("M")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),Function(self.domain),1,1))
    def test_setCoefficient_A_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numarray.ones((d,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,d),Function(self.domain),1,1))
    def test_setCoefficient_B_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numarray.ones((d,)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_C_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numarray.ones((d,)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),Function(self.domain),1,1))
    def test_setCoefficient_D_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=1.)
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),Function(self.domain),1,1))
    def test_setCoefficient_X_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numarray.ones((d,)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((d,),Function(self.domain),1))
    def test_setCoefficient_Y_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=1.)
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),Function(self.domain),1))
    def test_setCoefficient_y_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=1.)
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1))
    def test_setCoefficient_d_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=1.)
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_m_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=1.)
        coeff=mypde.getCoefficient("m")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),FunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=1.)
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=1.)
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),FunctionOnContactZero(self.domain),1))

    def test_setCoefficient_M_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M_reduced=1.)
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numarray.ones((d,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numarray.ones((d,)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numarray.ones((d,)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=1.)
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numarray.ones((d,)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=1.)
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=1.)
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_m_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m_reduced=1.)
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=1.)
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=1.)
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=1.)
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    def test_setCoefficient_r_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_q_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),Solution(self.domain),1))
    def test_setCoefficient_r_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=1.)
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))
    def test_setCoefficient_q_Scalar_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=1.)
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((),ReducedSolution(self.domain),1))

    def test_setCoefficient_M_reduced_Scalar_usingM(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_A_reduced_Scalar_usingA(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numarray.ones((d,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,d),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_B_reduced_Scalar_usingB(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numarray.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_C_reduced_Scalar_usingC(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numarray.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_D_reduced_Scalar_usingD(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1,1))
    def test_setCoefficient_X_reduced_Scalar_usingX(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=Data(numarray.ones((d,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((d,),ReducedFunction(self.domain),1))
    def test_setCoefficient_Y_reduced_Scalar_usingY(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Scalar(1.,ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunction(self.domain),1))
    def test_setCoefficient_y_reduced_Scalar_using_y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))
    def test_setCoefficient_m_reduced_Scalar_using_m(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_reduced_Scalar_using_d(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1,1))
    def test_setCoefficient_d_contact_reduced_Scalar_using_d_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1,1))
    def test_setCoefficient_y_contact_reduced_Scalar_using_y_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Scalar(1.,ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnContactZero(self.domain),1))
    #
    #  set coefficients for systems:
    #
    def test_setCoefficient_M_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("M")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_A_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=numarray.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_B_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=numarray.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_C_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=numarray.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N,d),Function(self.domain),self.N,self.N))
    def test_setCoefficient_D_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),Function(self.domain),self.N,self.N))
    def test_setCoefficient_X_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=numarray.ones((self.N,d)))
        coeff=mypde.getCoefficient("X")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,d),Function(self.domain),self.N))
    def test_setCoefficient_Y_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("Y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),Function(self.domain),self.N))
    def test_setCoefficient_y_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),FunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_m_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("m")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),FunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),FunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),FunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_M_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A_reduced=numarray.ones((self.N,d,self.N,d)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B_reduced=numarray.ones((self.N,d,self.N)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C_reduced=numarray.ones((self.N,self.N,d)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X_reduced=numarray.ones((self.N,d)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y_reduced=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_System_reduced(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_reduced=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_m_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact_reduced=numarray.ones((self.N,self.N)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact_reduced=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))
    def test_setCoefficient_r_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(r=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_q_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setValue(q=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),Solution(self.domain),self.N))
    def test_setCoefficient_r_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(r=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("r")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))
    def test_setCoefficient_q_System_reducedOn(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numEquations=3,debug=self.DEBUG)
        mypde.setReducedOrderOn()
        mypde.setValue(q=numarray.ones((self.N,)))
        coeff=mypde.getCoefficient("q")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions()),((self.N,),ReducedSolution(self.domain),self.N))

    def test_setCoefficient_M_System_reduced_using_D(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(M=Data(numarray.ones((self.N,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("M_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_A_reduced_System_using_A(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(A=Data(numarray.ones((self.N,d,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("A_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_B_reduced_System_using_B(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(B=Data(numarray.ones((self.N,d,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("B_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,d,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_C_reduced_System_using_C(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(C=Data(numarray.ones((self.N,self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("C_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N,d),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_D_System_reduced_using_D(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(D=Data(numarray.ones((self.N,self.N)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("D_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunction(self.domain),self.N,self.N))
    def test_setCoefficient_X_System_reduced_using_X(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(X=Data(numarray.ones((self.N,d)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("X_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,d),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_Y_System_reduced_using_Y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(Y=Data(numarray.ones((self.N,)),ReducedFunction(self.domain)))
        coeff=mypde.getCoefficient("Y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunction(self.domain),self.N))
    def test_setCoefficient_y_reduced_System_using_y(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Data(numarray.ones((self.N,)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnBoundary(self.domain),self.N))
    def test_setCoefficient_m_reduced_System_using_m(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(m=Data(numarray.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("m_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_reduced_System_using_d(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d=Data(numarray.ones((self.N,self.N)),ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficient("d_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnBoundary(self.domain),self.N,self.N))
    def test_setCoefficient_d_contact_reduced_System_using_d_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        mypde.setValue(d_contact=Data(numarray.ones((self.N,self.N)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("d_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumSolutions(),mypde.getNumEquations()),((self.N,self.N),ReducedFunctionOnContactZero(self.domain),self.N,self.N))
    def test_setCoefficient_y_contact_reduced_System_using_y_contact(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y_contact=Data(numarray.ones((self.N,)),ReducedFunctionOnContactZero(self.domain)))
        coeff=mypde.getCoefficient("y_contact_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((self.N,),ReducedFunctionOnContactZero(self.domain),self.N))

    def test_symmetryCheckTrue_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=100*numarray.ones((self.N,self.N))
        A=numarray.ones((self.N,d,self.N,d))
        C=2*numarray.ones((self.N,self.N,d))
        B=2*numarray.ones((self.N,d,self.N))
        D=3*numarray.ones((self.N,self.N))
        d=4*numarray.ones((self.N,self.N))
        m=64*numarray.ones((self.N,self.N))
        d_contact=5*numarray.ones((self.N,self.N))
        mypde.setValue(M=M,A=A,B=B,C=C,D=D,d=d,d_contact=d_contact,m=m,M_reduced=-M,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact, m_reduced=-m)
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_M_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=numarray.ones((self.N,self.N))
        M[1,0]=0.
        mypde.setValue(M=M)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_BC_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((self.N,self.N,d))
        B=2*numarray.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        D=3*numarray.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_m_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        m=4*numarray.ones((self.N,self.N))
        m[0,1]=0.
        mypde.setValue(m=m)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d=4*numarray.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numarray.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_M_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=3*numarray.ones((self.N,self.N))
        M[0,1]=0.
        mypde.setValue(M_reduced=M)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((self.N,d,self.N,d))
        A[1,1,1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_BC_reduced_System(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((self.N,self.N,d))
        B=2*numarray.ones((self.N,d,self.N))
        B[0,0,1]=1.
        mypde.setValue(B_reduced=B,C_reduced=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_D_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        D=3*numarray.ones((self.N,self.N))
        D[0,1]=0.
        mypde.setValue(D_reduced=D)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_m_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        m=4*numarray.ones((self.N,self.N))
        m[0,1]=0.
        mypde.setValue(m_reduced=m)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d=4*numarray.ones((self.N,self.N))
        d[0,1]=0.
        mypde.setValue(d_reduced=d)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_d_contact_reduced_System(self):
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        d_contact=5*numarray.ones((self.N,self.N))
        d_contact[0,1]=0.
        mypde.setValue(d_contact_reduced=d_contact)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")

    #==============================================================
    def test_symmetryCheckTrue_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        M=100
        A=numarray.ones((d,d))
        C=2*numarray.ones((d,))
        B=2*numarray.ones((d,))
        D=3
        m=10
        d=4
        d_contact=5
        mypde.setValue(M=M,A=A,B=B,C=C,D=D,d=d,m=m,d_contact=d_contact,M_reduced=-M,A_reduced=-A,B_reduced=-B,C_reduced=-C,D_reduced=-D,d_reduced=-d,d_contact_reduced=-d_contact,m_reduced=-m)
        self.failUnless(mypde.checkSymmetry(verbose=False),"symmetry detected")

    def test_symmetryCheckFalse_A_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((d,))
        B=2*numarray.ones((d,))
        B[0]=1.
        mypde.setValue(B=B,C=C)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_A_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        A=numarray.ones((d,d))
        A[1,0]=0.
        mypde.setValue(A_reduced=A)
        self.failUnless(not mypde.checkSymmetry(verbose=False),"symmetry detected")
    def test_symmetryCheckFalse_BC_reduced_Scalar(self):
        d=self.domain.getDim()
        mypde=TransportPDE(self.domain,debug=self.DEBUG)
        C=2*numarray.ones((d,))
        B=2*numarray.ones((d,))
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
