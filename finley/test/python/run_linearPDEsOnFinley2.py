
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

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_linearPDEs import Test_Poisson,Test_LinearPDE, Test_LinearPDE_noLumping, Test_TransportPDE
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_2Do2, Test_assemblage_3Do1, Test_assemblage_3Do2, \
                            Test_assemblage_2Do1_Contact,Test_assemblage_2Do2_Contact, Test_assemblage_3Do1_Contact, Test_assemblage_3Do2_Contact
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import *
from esys.finley import Rectangle,Brick,JoinFaces, ReadMesh
import sys


try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

NE=6 # number of element in each spatial direction (must be even)

class Test_LinearPDEOnFinleyTet2DOrder1(Test_LinearPDE,Test_pdetools,Test_assemblage_2Do1, Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet2DOrder2(Test_LinearPDE_noLumping,Test_pdetools_noLumping,Test_assemblage_2Do2, Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order2.fly"),optimize=False)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet3DOrder1(Test_LinearPDE,Test_pdetools,Test_assemblage_3Do1, Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet3DOrder2(Test_LinearPDE,Test_pdetools,Test_assemblage_3Do2, Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order2.fly"),optimize=False)
        self.order = 2
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

