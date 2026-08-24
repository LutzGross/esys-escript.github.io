
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

from test_simplesolve import ComplexSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_TRILINOS = hasFeature('trilinos')
skip_muelu_long = False #hasFeature("longindex")

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
mpiSize=getMPISizeWorld()


@unittest.skipIf(not HAVE_TRILINOS, "Trilinos not available")
class ComplexSolveOnTrilinos(ComplexSolveTestCase):
    pass


## direct
# TODO
# class Test_ComplexSolveOxley2D_Trilinos_Direct(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Rectangle(n0=NE0, n1=NE1)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

### BiCGStab + Jacobi

# TODO
# class Test_ComplexSolveOxley2D_Trilinos_BICGSTAB_Jacobi(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Rectangle(n0=NE0, n1=NE1)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

### GMRES + Jacobi

#TODO
# class Test_ComplexSolveOxley2D_Trilinos_GMRES_Jacobi(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Rectangle(n0=NE0, n1=NE1)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

### PCG + Jacobi

# class Test_ComplexSolveOxley2D_Trilinos_PCG_Jacobi(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Rectangle(n0=NE0, n1=NE1)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# ### PCG + AMG

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveOxley2D_Trilinos_PCG_AMG(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Rectangle(n0=NE0, n1=NE1)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

### PCG + ILUT

# class Test_ComplexSolveOxley2D_Trilinos_PCG_ILUT(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Rectangle(n0=NE0, n1=NE1)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain


# # BRICK
# class Test_ComplexSolveOxley3D_Trilinos_Direct(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveOxley3D_Trilinos_BICGSTAB_Jacobi(ComplexSolveOnTrilinos):
#     SOLVER_TOL = 1.e-9
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveOxley3D_Trilinos_GMRES_Jacobi(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveOxley3D_Trilinos_PCG_Jacobi(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveOxley3D_Trilinos_PCG_AMG(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 3)

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveOxley3D_Trilinos_PCG_ILUT(ComplexSolveOnTrilinos):
#     SOLVER_TOL = 1.e-9
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain


# This suite currently runs nothing. It is kept, rather than deleted,
# because the disabled tests describe what the feature should do - but
# a suite that collects no tests reports success, so state the fact
# explicitly here and let it show up as a skip.
@unittest.skip("every Test_ComplexSolveOxley*_Trilinos_* subclass in this file is commented out, so complex Trilinos solves on oxley have no coverage.")
class Test_ComplexSolveOxley_Trilinos_Disabled(unittest.TestCase):
    def test_disabled(self):
        pass


if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)

