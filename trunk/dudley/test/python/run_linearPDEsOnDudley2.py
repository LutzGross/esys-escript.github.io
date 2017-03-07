
########################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
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

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_linearPDEs import Test_LinearPDE, Test_TransportPDE
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_3Do1
from test_pdetools import Test_pdetools
from esys.dudley import ReadMesh

try:
     DUDLEY_TEST_DATA=os.environ['DUDLEY_TEST_DATA']
except KeyError:
     DUDLEY_TEST_DATA='.'

DUDLEY_TEST_MESH_PATH=os.path.join(DUDLEY_TEST_DATA,"data_meshes")

class Test_LinearPDEOnDudleyTet2D(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnDudleyTet3D(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnDudleyTet2D(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnDudleyTet3D(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnDudleyTet2D(Test_assemblage_2Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnDudleyTet3D(Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnDudleyTet2D(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnDudleyTet3D(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

