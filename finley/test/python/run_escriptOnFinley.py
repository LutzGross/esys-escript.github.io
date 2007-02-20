# $Id:$

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
from test_objects import Test_Dump as Test_Dump

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'


NE=4 # number elements, must be even
class Test_DumpOnFinley(Test_Dump):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
       self.domain_with_different_number_of_samples =Rectangle(2*NE,NE+1,2)
       self.domain_with_different_number_of_data_points_per_sample =Rectangle(2*NE,NE+1,2,integrationOrder=2)
       self.domain_with_different_sample_ordering =Rectangle(1,(NE+1)*NE,2)
       self.filebase=FINLEY_WORKDIR

   def tearDown(self):
       del self.domain
       del self.domain_with_different_number_of_samples
       del self.domain_with_different_number_of_data_points_per_sample
       del self.domain_with_different_sample_ordering

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_DumpOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
   
