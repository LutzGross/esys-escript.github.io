
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for PDE solvers on dudley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.dudley import Rectangle, Brick
from esys.escript import hasFeature
from esys.escript.linearPDEs import SolverOptions

HAVE_PASO = hasFeature('paso')

# number of elements in the spatial directions
NE0=12
NE1=13
NE2=8
OPTIMIZE=True

@unittest.skipIf(not HAVE_PASO, "PASO not available")
class SimpleSolveOnPaso(SimpleSolveTestCase):
    pass

### BiCGStab + Jacobi

class Test_SimpleSolveDudleyRect_Paso_BICGSTAB_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_BICGSTAB_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### PCG + Jacobi

class Test_SimpleSolveDudleyRect_Paso_PCG_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_PCG_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### TFQMR + Jacobi

class Test_SimpleSolveDudleyRect_Paso_TFQMR_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_TFQMR_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### MINRES + Jacobi

class Test_SimpleSolveDudleyRect_Paso_MINRES_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_MINRES_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### BiCGStab + Gauss-Seidel

class Test_SimpleSolveDudleyRect_Paso_BICGSTAB_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_BICGSTAB_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### PCG + Gauss-Seidel

class Test_SimpleSolveDudleyRect_Paso_PCG_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_PCG_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### TFQMR + Gauss-Seidel

class Test_SimpleSolveDudleyRect_Paso_TFQMR_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_TFQMR_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### MINRES + Gauss-Seidel

class Test_SimpleSolveDudleyRect_Paso_MINRES_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_MINRES_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### BiCGStab + RILU

class Test_SimpleSolveDudleyRect_Paso_BICGSTAB_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_BICGSTAB_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### PCG + RILU

class Test_SimpleSolveDudleyRect_Paso_PCG_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_PCG_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### TFQMR + RILU

class Test_SimpleSolveDudleyRect_Paso_TFQMR_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_TFQMR_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### MINRES + RILU

class Test_SimpleSolveDudleyRect_Paso_MINRES_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveDudleyBrick_Paso_MINRES_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
