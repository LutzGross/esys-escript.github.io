# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import unittest
from test_util import Test_util as Test_util
from test_util import Test_Util_SpatialFunctions
from test_symbols import Test_symbols

from esys.escript import FunctionOnBoundary
from esys.finley import Rectangle,Brick,JoinFaces
import sys

class Test_UtilOnFinley(Test_util,Test_symbols):
   def setUp(self):
       self.domain =Rectangle(5,5,2)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain

class Test_Util_SpatialFunctionsOnFinley2DOrder1(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Rectangle(n0=3,n1=12,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Rectangle(n0=3,n1=12,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2])

class Test_Util_SpatialFunctionsOnFinley2DOrder2(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Rectangle(n0=3,n1=6,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Rectangle(n0=3,n1=6,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2])

class Test_Util_SpatialFunctionsOnFinley3DOrder1(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Brick(n0=6,n1=12,n2=12,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Brick(n0=6,n1=12,n2=12,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2])

class Test_Util_SpatialFunctionsOnFinley3DOrder2(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Brick(n0=3,n1=6,n2=6,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Brick(n0=3,n1=6,n2=6,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2])

if __name__ == '__main__':
   suite = unittest.TestSuite()
#   suite.addTest(unittest.makeSuite(Test_UtilOnFinley))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinley2DOrder1))
#   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinley2DOrder2))
#   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinley3DOrder1))
#   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinley3DOrder2))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
   
