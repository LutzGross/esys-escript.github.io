
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
from test_pdetools import Test_pdetools
from esys.escript import hasFeature
from esys.finley import Brick

NE=10 # number of element in each spatial direction (must be even)

# trilinos and paso are different
TOL = 1.e-7 if hasFeature('paso') else 5.e-7

class Test_PDEToolsOnFinleyHex3DOrder2(Test_pdetools):
   RES_TOL=TOL
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,2)
        self.order = 2
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

