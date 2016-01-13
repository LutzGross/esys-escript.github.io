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
from esys.finley import Rectangle, JoinFaces, Brick

import numarray
FINLEY_TEST_MESH_PATH="data_meshes/"

NE=6 # number of element in each spatial direction (must be even)

class Test_AssemblePDEwithFinley_3Do2_Contact_withElementsOnFace(unittest.TestCase):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.,0.])
       self.domain = JoinFaces([d1,d2])
       self.domain.write("m.fly")
   def test_assemblage_3D_solO2_coeffOFull_NEqu1_B_Const_typeContact_comp0(self):
    x=self.domain.getX()
    u=(-4)+3*x[2]+6*x[2]**2-2*x[1]-5*x[1]*x[2]-8*x[1]**2+3*x[0]+7*x[0]*x[2]-9*x[0]*x[1]-4*x[0]**2
    B_test=Data(0.,(3,),Function(self.domain))
    B_test[0]=1
    Y_test=(-3)-7*x[2]+9*x[1]+8*x[0]
    x_boundary=FunctionOnBoundary(self.domain).getX()
    n=whereZero(x_boundary[0]   ,self.ABS_TOL)*numarray.array([-1., 0., 0.])+whereZero(x_boundary[0]-1.,self.ABS_TOL)*numarray.array([ 1., 0., 0.])+whereZero(x_boundary[1]   ,self.ABS_TOL)*numarray.array([ 0.,-1., 0.])+whereZero(x_boundary[1]-1.,self.ABS_TOL)*numarray.array([ 0., 1., 0.])+whereZero(x_boundary[2]   ,self.ABS_TOL)*numarray.array([ 0., 0.,-1.])+whereZero(x_boundary[2]-1.,self.ABS_TOL)*numarray.array([ 0., 0., 1.])
    y_test=n[0]*((-4)+3*x[2]+6*x[2]**2-2*x[1]-5*x[1]*x[2]-8*x[1]**2+3*x[0]+7*x[0]*x[2]-9*x[0]*x[1]-4*x[0]**2)
    n_contact=FunctionOnContactZero(self.domain).getNormal()
    y_contact_test=n_contact[0]*(4-3*x[2]-6*x[2]**2+2*x[1]+5*x[1]*x[2]+8*x[1]**2-3*x[0]-7*x[0]*x[2]+9*x[0]*x[1]+4*x[0]**2)
    pde=LinearPDE(self.domain)
    pde.setValue(B=B_test, Y=Y_test, y=y_test, y_contact=y_contact_test)
    r=pde.getResidual(u)
    rhs=pde.getRightHandSide()
    self.failUnless(Lsup(rhs)>0,"right hand side is zero")
    self.failUnless(Lsup(r)<=self.RES_TOL*Lsup(rhs),"residual is too big")
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_3Do2_Contact_withElementsOnFace))
   s=unittest.TextTestRunner(verbosity=2).run(suite)


