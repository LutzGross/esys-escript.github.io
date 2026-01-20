
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test suite for PDE solvers on finley
"""

from test_simplesolve import ComplexSolveTestCase, ComplexSolveTestCaseOrder2
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import Data, Solution, Vector, hasFeature, getMPISizeWorld
from esys.finley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_MUMPS = hasFeature('mumps')
skip_amg = True
skip_muelu_long = False #hasFeature("longindex")
mpiSize = getMPISizeWorld()

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
OPTIMIZE=True

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(mpiSize > 1, "MUMPS support for single MPI rank only")
class ComplexSolveOnMumps(ComplexSolveTestCase):
    pass
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(mpiSize > 1, "MUMPS support for single MPI rank only")
class ComplexSolveOnMumpsOrder2(ComplexSolveTestCaseOrder2):
    pass

## direct
class Test_ComplexSolveFinleyRect_Order1_Mumps_Direct(ComplexSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_ComplexSolveFinleyRect_Order2_Mumps_Direct(ComplexSolveOnMumpsOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_ComplexSolveFinleyBrick_Order1_Mumps_Direct(ComplexSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_ComplexSolveFinleyBrick_Order2_Mumps_Direct(ComplexSolveOnMumpsOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
        

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

