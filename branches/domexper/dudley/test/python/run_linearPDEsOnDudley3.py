
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

"""
Test suite for the linearPDE  and pdetools test on dudley

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
from test_linearPDEs import Test_Poisson,Test_LinearPDE, Test_LinearPDE_noLumping, Test_TransportPDE, Test_Helmholtz, Test_LameEquation
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_2Do2, Test_assemblage_3Do1, Test_assemblage_3Do2, \
                            Test_assemblage_2Do1_Contact,Test_assemblage_2Do2_Contact, Test_assemblage_3Do1_Contact, Test_assemblage_3Do2_Contact
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import *
from esys.dudley import Rectangle,Brick,JoinFaces, ReadMesh
import sys


try:
     DUDLEY_TEST_DATA=os.environ['DUDLEY_TEST_DATA']
except KeyError:
     DUDLEY_TEST_DATA='.'

DUDLEY_TEST_MESH_PATH=os.path.join(DUDLEY_TEST_DATA,"data_meshes")

NE=8 # number of element in each spatial direction (must be even)

class Test_LameOnDudley(Test_LameEquation):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2,useFullElementOrder=True)
   def tearDown(self):
        del self.domain
class Test_PoissonOnDudley(Test_Poisson):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2,useFullElementOrder=True)
   def tearDown(self):
        del self.domain

class Test_HelmholtzOnDudley(Test_Helmholtz):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2,useFullElementOrder=True)
   def tearDown(self):
        del self.domain



if __name__ == '__main__':
   suite = unittest.TestSuite()
   if True :
      suite.addTest(unittest.makeSuite(Test_PoissonOnDudley))
      suite.addTest(unittest.makeSuite(Test_HelmholtzOnDudley))
      suite.addTest(unittest.makeSuite(Test_LameOnDudley))
   else:
      pass

   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

