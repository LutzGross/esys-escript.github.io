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
   RES_TOL=1.e-2

   def setUp(self):
        self.order=1
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_2D_order2.fly")
   def tearDown(self):
        del self.order
        del self.domain
   def testProjector_rank0_fast(self):
      x=ContinuousFunction(self.domain).getX()
      h=Lsup(self.domain.getSize())
      p=Projector(self.domain,reduce=False,fast=True)
      td_ref=x[0]
      td=p(td_ref.interpolate(Function(self.domain)))
      print td
      print td_ref-td
      self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_util2))
   s=unittest.TextTestRunner(verbosity=2).run(suite)


