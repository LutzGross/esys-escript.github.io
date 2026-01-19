
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
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
class SimpleSolveOnMumps(SimpleSolveTestCase):
    pass
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(mpiSize > 1, "MUMPS support for single MPI rank only")
class SimpleSolveOnMumpsOrder2(SimpleSolveTestCaseOrder2):
    pass

## direct
class Test_SimpleSolveFinleyRect_Order1_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_SimpleSolveFinleyRect_Order2_Mumps_Direct(SimpleSolveOnMumpsOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_SimpleSolveFinleyBrick_Order1_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_SimpleSolveFinleyBrick_Order2_Mumps_Direct(SimpleSolveOnMumpsOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
        

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

