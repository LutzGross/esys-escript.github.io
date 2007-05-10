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
from esys.escript.pdetools import Projector
from esys.finley import Rectangle
# from test_pdetools import Test_pdetools

NE=6
import numarray
class Test_LinearPDEOnFinleyHex2DOrder1(unittest.TestCase):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,1)
   def tearDown(self):
        del self.domain
   def testProjector_rank3_fast_reduced_with_reduced_input(self):
      for i in range(800):
	      print i,i
	      f=ContinuousFunction(self.domain)
	      x=f.getX()
	      h=Lsup(self.domain.getSize())
	      p=Projector(self.domain, reduce=True,fast=True)
	      # td_ref=numarray.array([[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]])
      	      a=Data([1.,1.],ReducedFunction(self.domain))
              # Lsup(a)
	      td=p(a)
              # self.failUnless(Lsup(td-td_ref)<Lsup(td_ref)*h,"value wrong")

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyHex2DOrder1))
   s=unittest.TextTestRunner(verbosity=2).run(suite)


