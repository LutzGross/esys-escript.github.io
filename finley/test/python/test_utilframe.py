# $Id:$
"""
frame to ran a single test out of the Test_util suite
"""

import unittest
from esys.escript import *
from esys.finley import Rectangle
import numarray

class Test_util2(unittest.TestCase):
   RES_TOL=1.e-7
   def setUp(self):
       self.__dom =Rectangle(10,10,2)
       self.functionspace = FunctionOnBoundary(self.__dom) # due to a bug in escript python needs to hold a reference to the domain

   def test_wherePositive_array_rank0(self):
      arg=numarray.array(-15.0739210922)
      res=wherePositive(arg)
      ref=numarray.array(0.0)
      print res
      if isinstance(res,numarray.NumArray):
         self.failUnlessEqual(res.shape,(),"wrong shape of result.")
      else:
         self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")


if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_util2))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
