
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
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test suite for PDE solvers on oxley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_PASO = hasFeature('paso')

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
mpiSize=getMPISizeWorld()


@unittest.skipIf(not HAVE_PASO, "PASO not available")
class SimpleSolveOnPaso(SimpleSolveTestCase):
    pass


class Test_SimpleSolveOxley2D_Paso_DIRECT(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(n0=NE0, n1=NE1)
        self.method = SolverOptions.DIRECT
        # self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

# class Test_SimpleSolveOxley3D_Paso_DIRECT(SimpleSolveOnPaso):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain


class Test_SimpleSolveOxley2D_Paso_BICGSTAB_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(n0=NE0, n1=NE1)
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

# class Test_SimpleSolveOxley3D_Paso_BICGSTAB_Jacobi(SimpleSolveOnPaso):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

class Test_SimpleSolveOxley2D_Paso_PCG_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(n0=NE0, n1=NE1)
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

# class Test_SimpleSolveOxley3D_Paso_PCG_Jacobi(SimpleSolveOnPaso):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

class Test_SimpleSolveOxley2D_Paso_MINRES_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(n0=NE0, n1=NE1)
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

# class Test_SimpleSolveOxley3D_Paso_MINRES_Jacobi(SimpleSolveOnPaso):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.method = SolverOptions.MINRES
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

class Test_SimpleSolveOxley2D_Paso_TFQMR_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(n0=NE0, n1=NE1)
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

# class Test_SimpleSolveOxley3D_Paso_TFQMR_RILU(SimpleSolveOnPaso):
#     def setUp(self):
#         self.domain = Brick(n0=NE0, n1=NE1, n2=NE2)
#         self.method = SolverOptions.TFQMR
#         self.preconditioner = SolverOptions.RILU

#     def tearDown(self):
#         del self.domain


if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)

