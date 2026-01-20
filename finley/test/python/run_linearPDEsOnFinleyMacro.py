
########################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# Earth Systems Science Computational Center (ESSCC)
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
Earth Systems Science Computational Center (ESSCC)
https://github.com/LutzGross/esys-escript.github.io
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
from test_linearPDEs import Test_LinearPDE, Test_TransportPDE
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_3Do1
from test_pdetools import Test_pdetools
from esys.finley import Rectangle, Brick, ReadMesh

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

NE=3

class Test_AssemblageOnFinleyHex2DMacro(Test_assemblage_2Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyHex2DMacro(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyHex2DMacro(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyHex2DMacro(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnFinleyHex3DMacro(Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyHex3DMacro(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyHex3DMacro(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyHex3DMacro(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,-1,useElementsOnFace=0)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnFinleyTet2DMacro(Test_assemblage_2Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet2DMacro(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyTet2DMacro(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyTet2DMacro(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnFinleyTet3DMacro(Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet3DMacro(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyTet3DMacro(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyTet3DMacro(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_macro.fly"),optimize=False)
        self.order = 1
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

