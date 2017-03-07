
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
from esys.escript import *
from esys.dudley import Rectangle, Brick

NE=10 # number of element in each spatial direction (must be even)

class Test_LinearPDEOnDudleyRect(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnDudleyBrick(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnDudleyRect(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnDudleyBrick(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnDudleyRect(Test_assemblage_2Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_AssemblageOnDudleyBrick(Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnDudleyRect(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnDudleyBrick(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE)
        self.order = 1
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

