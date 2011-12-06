
########################################################
#
# Copyright (c) 2003-2011 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2011 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
from test_util import Test_util_no_tagged_data
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_symbols import Test_symbols

from esys.escript import *
from esys.ripley import Rectangle, Brick
import sys
import os

try:
     RIPLEY_TEST_DATA=os.environ['RIPLEY_TEST_DATA']
except KeyError:
     RIPLEY_TEST_DATA='.'

NE=4 # number elements

class Test_UtilOnRipley(Test_util_no_tagged_data,Test_symbols):
   def setUp(self):
       mpiSize=getMPISizeWorld()
       for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
           NX=x
           NY=mpiSize/x
           if NX*NY == mpiSize:
               break
       self.domain=Rectangle(n0=NE*NX, n1=NE*NY, d0=NX, d1=NY)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_Util_SpatialFunctionsOnRipleyRect(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        mpiSize=getMPISizeWorld()
        for x in [int(sqrt(mpiSize)),7,5,3,2,1]:
            NX=x
            NY=mpiSize/x
            if NX*NY == mpiSize:
                break
        self.domain = Rectangle(n0=NE*NX, n1=NE*NY, d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnRipleyBrick(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        mpiSize=getMPISizeWorld()
        for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
            NX=x[0]
            NY=x[1]
            NZ=mpiSize/(x[0]*x[1])
            if NX*NY*NZ == mpiSize:
                break
        self.domain = Brick(n0=NE*NX, n1=NE*NY, n2=NE*NZ, d0=NX, d1=NY, d2=NZ)
    def tearDown(self):
        del self.order
        del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   if True:
      suite.addTest(unittest.makeSuite(Test_UtilOnRipley))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnRipleyRect))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnRipleyBrick))
   else:
      pass
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

