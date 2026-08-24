
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
Test suite for PDE solvers on Oxley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

SOLVER="mkl"
HAVE_REQUESTED_SOLVER = hasFeature(SOLVER)

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
mpiSize=getMPISizeWorld()

    


@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "MKL runs on single rank only.")
class Test_SimpleSolveOxley2D_MKL(SimpleSolveTestCase):
    def setUp(self):
        self.domain = Rectangle(n0=NE0, n1=NE1)
        self.package = SolverOptions.MKL
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
        
@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "MKL runs on single rank only.")
class Test_SimpleSolveOxley3D_MKL(SimpleSolveTestCase):
    def setUp(self):
        self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
        self.package = SolverOptions.MKL
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain


if __name__ == "__main__":
    run_tests(__name__, exit_on_failure=True)
