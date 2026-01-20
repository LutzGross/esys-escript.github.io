
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test suite for PDE solvers on finley
"""

from test_simplesolve import SimpleSolveTestCase, SimpleSolveTestCaseOrder2
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, Data, Solution, Vector, hasFeature
from esys.finley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions



SOLVER='umfpack'
HAVE_REQUESTED_SOLVER = hasFeature(SOLVER) and hasFeature('paso')


mpiSize=getMPISizeWorld()
# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
OPTIMIZE=True



@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "UMFPACK runs on single rank only.")
class Test_SimpleSolveFinleyRect_Order1_PasoUMFPACK(SimpleSolveTestCase):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.UMFPACK
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "UMFPACK runs on single rank only.")
class Test_SimpleSolveFinleyRect_Order2_PasoUMFPACK(SimpleSolveTestCaseOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.UMFPACK
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "UMFPACK runs on single rank only.")
class Test_SimpleSolveFinleyBrick_Order1_PasoUMFPACK(SimpleSolveTestCase):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.UMFPACK
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "UMFPACK runs on single rank only.")
class Test_SimpleSolveFinleyBrick_Order2_PasoUMFPACK(SimpleSolveTestCaseOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.UMFPACK
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
