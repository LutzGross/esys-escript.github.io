
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
from test_linearPDEs import Test_Helmholtz, Test_LameEquation, Test_Poisson
from test_assemblage import Test_assemblage_2Do1_Contact, Test_assemblage_2Do2_Contact, Test_assemblage_3Do1_Contact, Test_assemblage_3Do2_Contact
from esys.finley import Rectangle, ReadMesh

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

NE=8 # number of element in each spatial direction (must be even)

class Test_LameOnFinley(Test_LameEquation):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE, NE, 2, useElementsOnFace=0, useFullElementOrder=True)
   def tearDown(self):
        del self.domain

class Test_PoissonOnFinley(Test_Poisson):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE, NE, 2, useElementsOnFace=0, useFullElementOrder=True)
   def tearDown(self):
        del self.domain

class Test_HelmholtzOnFinley(Test_Helmholtz):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE, NE, 2, useElementsOnFace=0, useFullElementOrder=True)
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_2Do1_Contact(Test_assemblage_2Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do1_Contact.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_2Do2_Contact(Test_assemblage_2Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do2_Contact.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do1_Contact(Test_assemblage_3Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do1_Contact.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do2_Contact(Test_assemblage_3Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do2_Contact.fly"))
   def tearDown(self):
        del self.domain


class Test_AssemblePDEwithFinley_2Do1_Contact_withElementsOnFace(Test_assemblage_2Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do1_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_2Do2_Contact_withElementsOnFace(Test_assemblage_2Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do2_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do1_Contact_withElementsOnFace(Test_assemblage_3Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do1_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do2_Contact_withElementsOnFace(Test_assemblage_3Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do2_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

