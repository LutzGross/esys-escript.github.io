# $Id$

import sys
import unittest

from esys.escript.test_linearPDEs import Test_Poisson,Test_LinearPDE
from esys.finley import Rectangle,Brick

class Test_LinearPDEOnFinley2DOrder1(Test_LinearPDE):
    def setUp(self):
        self.domain = Rectangle(20,20,1)

class Test_LinearPDEOnFinley2DOrder2(Test_LinearPDE):
    def setUp(self):
        self.domain = Rectangle(20,20,2)

class Test_LinearPDEOnFinley3DOrder1(Test_LinearPDE):
    def setUp(self):
        self.domain = Brick(20,10,10,1)

class Test_LinearPDEOnFinley3DOrder2(Test_LinearPDE):
    def setUp(self):
        self.domain = Brick(10,10,20,2)

class Test_PoissonOnFinley(Test_Poisson):
    def setUp(self):
        self.domain = Rectangle(20,10,2)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley2DOrder1))
#   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley2DOrder2))
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley3DOrder1))
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley3DOrder2))
   suite.addTest(unittest.makeSuite(Test_PoissonOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
