
########################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for the linearPDE  and pdetools test on finley

:remark:

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import os
import unittest
from test_nonLinearPDE import Test_nonLinearPDEs, Test_nlpde
from esys.escript import *
from esys.dudley import Rectangle,Brick
import sys


class Test_nonLinearPDE(Test_nlpde):
   def setUp(self):
        self.domain = Rectangle(l0=10.,l1=10.,n0=10, n1=10) 
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   if True :
      suite.addTest(unittest.makeSuite(Test_nonLinearPDE))
   else:
      suite.addTest(Test_LinearPDEOnFinleyHex2DOrder1("testProjector_rank1_fast_reduced"))
      pass

   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

