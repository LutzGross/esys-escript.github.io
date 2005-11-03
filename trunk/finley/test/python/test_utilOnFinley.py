# $Id$

import unittest
from esys.escript.test_util import Test_Util
from esys.escript import FunctionOnBoundary
from esys.finley import Rectangle
import sys

class Test_UtilOnFinley(Test_Util):
   def setUp(self):
       self.__dom =Rectangle(10,10,2)
       self.functionspace = FunctionOnBoundary(self.__dom) # due to a bug in escript python needs to hold a reference to the domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_UtilOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
   
