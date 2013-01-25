#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

import unittest
from test_util import Test_util_no_tagged_data as Test_util
from test_symbols import Test_symbols
from esys.escript import ContinuousFunction
from esys.bruce import Rectangle
import sys

class Test_UtilOnBruce(Test_util,Test_symbols):
   def setUp(self):
       self.domain =Rectangle(10,10)
       self.functionspace = ContinuousFunction(self.domain)

   def tearDown(self):
    del self.domain
    del self.functionspace


if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_UtilOnBruce))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)
