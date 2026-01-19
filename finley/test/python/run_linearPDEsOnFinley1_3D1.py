
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

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_linearPDEs import Test_LinearPDE, Test_TransportPDE
from test_assemblage import Test_assemblage_3Do1
from test_pdetools import Test_pdetools
from esys.finley import Brick

NE=10 # number of element in each spatial direction (must be even)

class Test_AssemblageOnFinleyHex3DOrder1(Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,1)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyHex3DOrder1(Test_LinearPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,1)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_PDEToolsOnFinleyHex3DOrder1(Test_pdetools):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,1)
        self.order = 1
   def tearDown(self):
        del self.domain

class Test_TransportPDEOnFinleyHex3DOrder1(Test_TransportPDE):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,1)
        self.order = 1
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

