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
from esys.finley import ReadMesh
from esys.escript.pdetools import Projector
import numarray
FINLEY_TEST_MESH_PATH="data_meshes/"

class Test_util2(unittest.TestCase):
   RES_TOL=1.e-8
   def setUp(self):
        self.order=1
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_3D_order1.fly")
   def tearDown(self):
        del self.order
        del self.domain

   def test_normal_FunctionOnBoundary(self):
     """
     test getNormal() on boundary

     assumptions: FunctionOnBoundary(self.domain) exists
     """
     dim=self.domain.getDim()
     f=FunctionOnBoundary(self.domain)
     x=f.getX()
     ref=Vector(0.,what=f)
     if dim==3:
         ref+=whereZero(x[0]-1.,tol=self.RES_TOL)*[1,0,0]
         ref+=whereZero(x[0],tol=self.RES_TOL)*[-1,0,0]
         ref+=whereZero(x[1]-1.,tol=self.RES_TOL)*[0,1,0]
         ref+=whereZero(x[1],tol=self.RES_TOL)*[0,-1,0]
         ref+=whereZero(x[2]-1.,tol=self.RES_TOL)*[0,0,1]
         ref+=whereZero(x[2],tol=self.RES_TOL)*[0,0,-1]
     else:
         ref+=whereZero(x[0]-1.,tol=self.RES_TOL)*[1,0]
         ref+=whereZero(x[0],tol=self.RES_TOL)*[-1,0]
         ref+=whereZero(x[1]-1.,tol=self.RES_TOL)*[0,1]
         ref+=whereZero(x[1],tol=self.RES_TOL)*[0,-1]

     res=f.getNormal()
     print length(ref-res)
     self.failUnlessEqual(res.getShape(),(dim,),"wrong shape of result.")
     self.failUnlessEqual(res.getFunctionSpace(),f,"wrong functionspace of result.")
     self.failUnless(Lsup(ref-res)<=self.RES_TOL,"wrong result")

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_util2))
   s=unittest.TextTestRunner(verbosity=2).run(suite)


