# -*- coding: utf-8 -*-

########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
from test_util import Test_util as Test_util
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary, Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_symbols import Test_symbols

from esys.escript import *
from esys.buckley import Buckley
import sys
import os

try:
     BUCKLEY_TEST_DATA=os.environ['BUCKLEY_TEST_DATA']
except KeyError:
     BUCKLEY_TEST_DATA='.'

BUCKLEY_TEST_MESH_PATH=os.path.join(BUCKLEY_TEST_DATA,"data_meshes")


NR=4 # Number of refinements

class Test_UtilOnBuckley(Test_util,Test_symbols):
   def setUp(self):
       self.domain =Buckley(1,1,1)
       self.domain.refineAll(NR)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain




class Test_Util_SpatialFunctionsOnBuckley(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.domain = Buckley(1,1,1)
        self.domain.refineAll(NR)
    def tearDown(self):
        del self.domain


if __name__ == '__main__':
   suite = unittest.TestSuite()
   if True:
      suite.addTest(unittest.makeSuite(Test_UtilOnBuckley))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnBuckley))
   else:
      pass
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

