
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for PDE solvers on Dudley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, Data, Solution, Vector, hasFeature
from esys.dudley import Rectangle, Brick
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
class Test_SimpleSolveDudleyRect_Umfpack_Direct(SimpleSolveTestCase):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.UMFPACK
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
        
@unittest.skipIf(not HAVE_REQUESTED_SOLVER, "%s not available"%SOLVER)
@unittest.skipIf(mpiSize > 1, "UMFPACK runs on single rank only.")
class Test_SimpleSolveDudleyBrick_Umfpack_Direct(SimpleSolveTestCase):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.UMFPACK
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
        

