
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

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-dudley"

"""
Test suite for PDE solvers on dudley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import Data, Solution, Vector, hasFeature
from esys.dudley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_MUMPS = hasFeature('mumps')
skip_amg = True
skip_muelu_long = False #hasFeature("longindex")

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
OPTIMIZE=True

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class SimpleSolveOnMumps(SimpleSolveTestCase):
    pass

## direct
class Test_SimpleSolveDudleyRect_Order1_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi

class Test_SimpleSolveDudleyRect_Order1_Mumps_BICGSTAB_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

## direct
class Test_SimpleSolveDudleyBrick_Order1_Mumps_Direct(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_SimpleSolveDudleyBrick_Order1_Mumps_BICGSTAB_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### PCG + Jacobi

class Test_SimpleSolveDudleyRect_Order1_Mumps_PCG_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_PCG_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### TFQMR + Jacobi

class Test_SimpleSolveDudleyRect_Order1_Mumps_TFQMR_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_TFQMR_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### MINRES + Jacobi

class Test_SimpleSolveDudleyRect_Order1_Mumps_MINRES_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_MINRES_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### BiCGStab + Gauss-Seidel

class Test_SimpleSolveDudleyRect_Order1_Mumps_BICGSTAB_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_BICGSTAB_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### PCG + AMG

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveDudleyRect_Order1_Mumps_PCG_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveDudleyBrick_Order1_Mumps_PCG_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

### TFQMR + Gauss-Seidel

class Test_SimpleSolveDudleyRect_Order1_Mumps_TFQMR_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_TFQMR_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### MINRES + AMG

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveDudleyRect_Order1_Mumps_MINRES_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveDudleyBrick_Order1_Mumps_MINRES_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

### BiCGStab + RILU

class Test_SimpleSolveDudleyRect_Order1_Mumps_BICGSTAB_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_BICGSTAB_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### PCG + RILU

class Test_SimpleSolveDudleyRect_Order1_Mumps_PCG_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_PCG_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### TFQMR + RILU

class Test_SimpleSolveDudleyRect_Order1_Mumps_TFQMR_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_TFQMR_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### MINRES + RILU

class Test_SimpleSolveDudleyRect_Order1_Mumps_MINRES_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_MINRES_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### PCG + ILUT

class Test_SimpleSolveDudleyRect_Order1_Mumps_PCG_ILUT(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Order1_Mumps_PCG_ILUT(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

