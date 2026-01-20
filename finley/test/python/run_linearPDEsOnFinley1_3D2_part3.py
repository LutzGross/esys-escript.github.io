
########################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test suite for the linearPDE and pdetools test on finley
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_assemblage import Test_assemblage_3Do2
from esys.finley import Brick

NE=10 # number of element in each spatial direction (must be even)

class Test_AssemblageOnFinleyHex3DOrder2(Test_assemblage_3Do2):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,2)
        self.order = 2
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

