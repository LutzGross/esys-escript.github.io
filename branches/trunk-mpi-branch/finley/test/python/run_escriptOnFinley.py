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
import tempfile

from esys.escript import *
from esys.finley import Rectangle
import sys
import os
from test_objects import Test_Dump, Test_SetDataPointValue
from test_objects import Test_Domain

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

NE=4 # number elements, must be even
class Test_DomainOnFinley(Test_Domain):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
   def tearDown(self):
       del self.domain
class Test_DataOpsOnFinley(Test_Dump): # , Test_SetDataPointValue):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
       self.domain_with_different_number_of_samples =Rectangle(2*NE,NE+1,2)
       self.domain_with_different_number_of_data_points_per_sample =Rectangle(2*NE,NE+1,2,integrationOrder=2)
       self.domain_with_different_sample_ordering =Rectangle(1,(NE+1)*NE,2)
       self.filename_base=FINLEY_WORKDIR

   def tearDown(self):
       del self.domain
       del self.domain_with_different_number_of_samples
       del self.domain_with_different_number_of_data_points_per_sample
       del self.domain_with_different_sample_ordering

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_DataOpsOnFinley))
   suite.addTest(unittest.makeSuite(Test_DomainOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
