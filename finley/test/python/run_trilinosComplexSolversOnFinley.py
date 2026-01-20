
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

from test_simplesolve import ComplexSolveTestCase, ComplexSolveTestCaseOrder2
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import Data, Solution, Vector, hasFeature
from esys.finley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_TRILINOS = hasFeature('trilinos')
skip_muelu_long = False #hasFeature("longindex")

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
OPTIMIZE=True

@unittest.skipIf(not HAVE_TRILINOS, "Trilinos not available")
class ComplexSolveOnTrilinos(ComplexSolveTestCase):
    pass
@unittest.skipIf(not HAVE_TRILINOS, "Trilinos not available")
class ComplexSolveOnTrilinosOrder2(ComplexSolveTestCaseOrder2):
    pass

## direct
class Test_ComplexSolveFinleyRect_Order1_Trilinos_Direct(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_ComplexSolveFinleyRect_Order2_Trilinos_Direct(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_ComplexSolveFinleyBrick_Order1_Trilinos_Direct(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

## direct
class Test_ComplexSolveFinleyBrick_Order2_Trilinos_Direct(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain
        
### BiCGStab + Jacobi
@unittest.skip("fails with Nan during iteration.")
class Test_ComplexSolveFinleyRect_Order1_Trilinos_BICGSTAB_Jacobi(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

@unittest.skip("fails with Nan during iteration.")
class Test_ComplexSolveFinleyRect_Order2_Trilinos_BICGSTAB_Jacobi(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain


        
class Test_ComplexSolveFinleyBrick_Order1_Trilinos_BICGSTAB_Jacobi(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain
        
@unittest.skip("convergence problems")
class Test_ComplexSolveFinleyBrick_Order2_Trilinos_BICGSTAB_Jacobi(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI
        self.REL_TOL=5.e-6
    def tearDown(self):
        del self.domain

### PCG + Jacobi

class Test_ComplexSolveFinleyRect_Order1_Trilinos_PCG_Jacobi(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyRect_Order2_Trilinos_PCG_Jacobi(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order1_Trilinos_PCG_Jacobi(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order2_Trilinos_PCG_Jacobi(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain
### BiCGStab + Gauss-Seidel

class Test_ComplexSolveFinleyRect_Order1_Trilinos_BICGSTAB_GaussSeidel(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

@unittest.skip("fails with Nan during iteration.")
class Test_ComplexSolveFinleyRect_Order2_Trilinos_BICGSTAB_GaussSeidel(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order1_Trilinos_BICGSTAB_GaussSeidel(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

@unittest.skip("convergence problems")
class Test_ComplexSolveFinleyBrick_Order2_Trilinos_BICGSTAB_GaussSeidel(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL
        self.REL_TOL=5.e-6

    def tearDown(self):
        del self.domain

### PCG + AMG

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveFinleyRect_Order1_Trilinos_PCG_AMG(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveFinleyRect_Order2_Trilinos_PCG_AMG(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveFinleyBrick_Order1_Trilinos_PCG_AMG(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveFinleyBrick_Order2_Trilinos_PCG_AMG(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain


### BiCGStab + RILU

class Test_ComplexSolveFinleyRect_Order1_Trilinos_BICGSTAB_RILU(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyRect_Order2_Trilinos_BICGSTAB_RILU(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order1_Trilinos_BICGSTAB_RILU(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order2_Trilinos_BICGSTAB_RILU(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### PCG + RILU

class Test_ComplexSolveFinleyRect_Order1_Trilinos_PCG_RILU(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyRect_Order2_Trilinos_PCG_RILU(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order1_Trilinos_PCG_RILU(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order2_Trilinos_PCG_RILU(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### PCG + ILUT

class Test_ComplexSolveFinleyRect_Order1_Trilinos_PCG_ILUT(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyRect_Order2_Trilinos_PCG_ILUT(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order1_Trilinos_PCG_ILUT(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveFinleyBrick_Order2_Trilinos_PCG_ILUT(ComplexSolveOnTrilinosOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

