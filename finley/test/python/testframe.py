# $Id:$
"""
frame to ran a single test out of the Test_util suite
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import unittest
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle, JoinFaces
import numarray
FINLEY_TEST_MESH_PATH="data_meshes/"

NE=6 # number of element in each spatial direction (must be even)
class Test_util2(unittest.TestCase):
  RES_TOL=1.e-7
  ABS_TOL=1.e-8
  def setUp(self):
       d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.])
       self.domain = JoinFaces([d1,d2])
  def tearDown(self):
        del self.domain

  #==================================================
  def test_assemblage_2D_solO1_coeffOFull_NEqu1_d_contact_Const_typeContact(self):
    x=self.domain.getX()
    jump=Data(0.,(),ContinuousFunction(self.domain))
    jump.setTaggedValue(2,1.)
    u=((-6)+x[1]-6*x[0])*jump
    d_contact_test=Data(4,(),FunctionOnContactZero(self.domain))
    y_contact_test=(-24)+4*x[1]-24*x[0]
    pde=LinearPDE(self.domain)
    pde.setValue(d_contact=d_contact_test, y_contact=y_contact_test)
    r=pde.getResidual(u)
    rhs=pde.getRightHandSide()
    self.failUnless(Lsup(rhs)>0,"right hand side is zero")
    self.failUnless(Lsup(r)<=self.RES_TOL*Lsup(rhs),"residual is too big")

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_util2))
   s=unittest.TextTestRunner(verbosity=2).run(suite)


