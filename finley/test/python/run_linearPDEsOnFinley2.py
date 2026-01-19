
########################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

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
from test_linearPDEs import Test_LinearPDE, Test_LinearPDE_noLumping, Test_TransportPDE
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_2Do2, Test_assemblage_3Do1, Test_assemblage_3Do2
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import hasFeature
from esys.finley import ReadMesh

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

# paso and trilinos are different
TOL = 1.e-7 if hasFeature('paso') else 5.e-7

class Test_AssemblageOnFinleyTet2DOrder1(Test_assemblage_2Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet2DOrder1(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyTet2DOrder1(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyTet2DOrder1(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnFinleyTet2DOrder2(Test_assemblage_2Do2):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order2.fly"),optimize=False)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet2DOrder2(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order2.fly"),optimize=False)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyTet2DOrder2(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order2.fly"),optimize=False)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyTet2DOrder2(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order2.fly"),optimize=False)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_AssemblageOnFinleyTet3DOrder1(Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order1.fly"), optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet3DOrder1(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyTet3DOrder1(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyTet3DOrder1(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=True)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnFinleyTet3DOrder2(Test_assemblage_3Do2):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order2.fly"),optimize=True)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet3DOrder2(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order2.fly"),optimize=True)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyTet3DOrder2(Test_pdetools):
   RES_TOL=TOL
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order2.fly"),optimize=True)
        self.order = 2
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyTet3DOrder2(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order2.fly"),optimize=True)
        self.order = 2
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

