
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

"""
Test suite for PDE solvers on Oxley multiresolution domains
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
# from esys.oxley import MultiResolutionDomain
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
for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
    NX=x
    NY=mpiSize//x
    if NX*NY == mpiSize:
        break

for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,2),(1,1)]:
    NXb=x[0]
    NYb=x[1]
    NZb=mpiSize//(x[0]*x[1])
    if NXb*NYb*NZb == mpiSize:
        break

def test_Rectangle(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    return m

# def Brick(**kwargs):
#     m = MultiResolutionDomain(3, **kwargs)
#     return m.getLevel(1)

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(mpiSize > 1, "MUMPS support for single MPI rank only")
class SimpleSolveOnMumps(SimpleSolveTestCase):
    pass

## MUMPS direct solver tests
## Note: MUMPS only supports direct solving - iterative methods are not available

@unittest.skip("Oxley matrix distribution issue - see issue #118")
class Test_SimpleSolveMultires2D_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

@unittest.skip("Oxley matrix distribution issue - see issue #118")
class Test_SimpleSolveMultires3D_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

