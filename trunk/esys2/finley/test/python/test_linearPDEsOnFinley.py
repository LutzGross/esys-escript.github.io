# $Id$

import unittest
from esys.escript.test_linearPDEs import Test_Poisson,Test_LinearPDE
from esys.finley import Rectangle
import sys

class Test_LinearPDEOnFinley(Test_LinearPDE):
    def setUp(self):
        self.domain = Rectangle(10,10,2)

class Test_PoissonOnFinley(Test_Poisson):
    def setUp(self):
        self.domain = Rectangle(10,10,2)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley))
   suite.addTest(unittest.makeSuite(Test_PoissonOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
   
