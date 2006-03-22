# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__licence__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licences/osl-3.0.php"""
import sys
import unittest
import os

from esys.escript import *
from esys import bruce
from esys.escript.test_symbols import Test_symbols
from esys.escript.test_util import Test_util

import numarray

class Test_symbols_OnBruce(Test_symbols):
    def setUp(self):
        self.domain = bruce.Brick(5,15,8)
        self.functionspace=escript.ContinuousFunction(self.domain)

class Test_util_OnBruce(Test_util):
    def setUp(self):
        self.domain = bruce.Brick(5,15,8)
        self.functionspace=escript.ContinuousFunction(self.domain)


class escriptTestCase(unittest.TestCase):

  def testPow(self):
    """Test the pow function."""
    myFuncSpac=escript.FunctionSpace()
    myVector=numarray.array([[1,2],[3,4]])
    myData=escript.Data(myVector,myFuncSpac,True)
    print myData**3

  def testFunctionSpace(self):
    """Test the creation of FunctionSpace objects."""
    mesh=bruce.Brick(1,1,1)
    cFunc=escript.ContinuousFunction(mesh)

  def testDataOperations(self):
    """Test the operations that may be performed on Data."""
    myFuncSpac=escript.FunctionSpace()
    myVector=numarray.array([[1,2],[3,4]])
    myData=escript.Data(myVector,myFuncSpac,True)
    myData2=3+myData
    print myData2
    myList=[[1,2],[3,4]]
    myData3=myList+myData
    print myData3
    myData3=myData+myData2+myData3
    print myData3
    myData4=myList-myData
    print myData4
    myData5=0-myData
    print myData5
    myData6=1/myData
    print myData6
    myData7=5*myData
    print myData7
    myData8=(3*myData)/3+2*myData-2*myData
    print myData8
    myData9=myData.sin()
    print myData9
    myData9=myData.cos()
    print myData9
    print myData9.wherePositive()
    print myData9.whereNegative()
    
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(escriptTestCase))
suite.addTest(unittest.makeSuite(Test_symbols_OnBruce))
suite.addTest(unittest.makeSuite(Test_util_OnBruce))
if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
  sys.exit(0)
else:
  sys.exit(1)
