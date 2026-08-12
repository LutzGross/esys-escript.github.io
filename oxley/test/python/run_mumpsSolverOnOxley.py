
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

"""
Test suite for PDE solvers on Oxley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_MUMPS = hasFeature('mumps')
skip_amg = True
skip_muelu_long = False #hasFeature("longindex")

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
OPTIMIZE=True
mpiSize=getMPISizeWorld()


@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(mpiSize > 1, "MUMPS support for single MPI rank only")
class SimpleSolveOnMumps(SimpleSolveTestCase):
    pass

## MUMPS direct solver tests
## Note: MUMPS only supports direct solving - iterative methods are not available

class Test_SimpleSolveOxley2D_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(n0=NE0, n1=NE1)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

class Test_SimpleSolveOxley3D_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

