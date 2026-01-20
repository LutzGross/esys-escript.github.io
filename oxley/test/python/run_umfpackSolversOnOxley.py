
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
Test suite for PDE solvers on Oxley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

SOLVER='umfpack'
HAVE_REQUESTED_SOLVER = hasFeature(SOLVER) and hasFeature('paso')

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
mpiSize=getMPISizeWorld()
for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
    NX=x
    NY=mpiSize//x
    if NX*NY == mpiSize:
        break

for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
    NXb=x[0]
    NYb=x[1]
    NZb=mpiSize//(x[0]*x[1])
    if NXb*NYb*NZb == mpiSize:
        break
    


@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "UMFPACK runs on single rank only.")
@unittest.skip("Oxley Rectangle meshes have SystemMatrixPattern errors with UMFPACK - see issue #118")
class Test_SimpleSolveOxley2D_UMFPACK(SimpleSolveTestCase):
    def setUp(self):
        self.domain = Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.UMFPACK
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
        
# @unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
# @unittest.skipIf(mpiSize > 1, "UMFPACK runs on single rank only.")
# class Test_SimpleSolveOxley3D_UMFPACK(SimpleSolveTestCase):
#     def setUp(self):
#         self.domain = Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
#         self.package = SolverOptions.UMFPACK
#         self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
