# $Id$

"""
Test suite for the linearPDE iand pdetools test on finley

@remark:

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Lutz Gross, l.gross@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"


import unittest
from esys.escript.test_linearPDEs import Test_Poisson,Test_LinearPDE
from esys.escript.test_pdetools import Test_pdetools
from esys.finley import Rectangle,Brick
import sys

class Test_LinearPDEOnFinley2DOrder1(Test_LinearPDE,Test_pdetools):
    def setUp(self):
        self.domain = Rectangle(50,50,1)

class Test_LinearPDEOnFinley2DOrder2(Test_LinearPDE,Test_pdetools):
    def setUp(self):
        self.domain = Rectangle(50,50,2)

class Test_LinearPDEOnFinley3DOrder1(Test_LinearPDE,Test_pdetools):
    def setUp(self):
        self.domain = Brick(20,10,10,1)

class Test_LinearPDEOnFinley3DOrder2(Test_LinearPDE,Test_pdetools):
    def setUp(self):
        self.domain = Brick(10,10,20,2)

class Test_PoissonOnFinley(Test_Poisson):
    def setUp(self):
        self.domain = Rectangle(20,10,2)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley2DOrder1))
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley2DOrder2))
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley3DOrder1))
   suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinley3DOrder2))
   suite.addTest(unittest.makeSuite(Test_PoissonOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
   
