
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
Test suite for PDE solvers on finley
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import Data, Solution, Vector, hasFeature
from esys.finley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_PASO = hasFeature('paso')

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
OPTIMIZE=True

@unittest.skipIf(not HAVE_PASO, "PASO not available")
class SimpleSolveOnPaso(SimpleSolveTestCase):
    pass

class SimpleSolveOrder2(SimpleSolveOnPaso):
    def _getGrad(self, system):
        """returns exact gradient"""
        dim = self.domain.getDim()
        x = Solution(self.domain).getX()
        if system:
            g_ex = Data(0., (dim,dim), Solution(self.domain))
            if dim == 2:
                g_ex[0,0] = 2.+8.*x[0]+ 5.*x[1]
                g_ex[0,1] = 3.+5.*x[0]+12.*x[1]
                g_ex[1,0] = 4.+2.*x[0]+ 6.*x[1]
                g_ex[1,1] = 2.+6.*x[0]+ 8.*x[1]
            else:
                g_ex[0,0] =  2.+6.*x[1]+8.*x[2]+18.*x[0]
                g_ex[0,1] =  3.+6.*x[0]+7.*x[2]+20.*x[1]
                g_ex[0,2] =  4.+7.*x[1]+8.*x[0]+22.*x[2]
                g_ex[1,0] =  4.+3.*x[1]-8.*x[2]- 4.*x[0]
                g_ex[1,1] =  1.+3.*x[0]+2.*x[2]+14.*x[1]
                g_ex[1,2] = -6.+2.*x[1]-8.*x[0]+10.*x[2]
                g_ex[2,0] =  7.-6.*x[1]+2.*x[2]+ 4.*x[0]
                g_ex[2,1] =  9.-6.*x[0]+8.*x[2]+16.*x[1]
                g_ex[2,2] =  2.+8.*x[1]+2.*x[0]+ 2.*x[2]
        else:
            g_ex = Data(0., (dim,), Solution(self.domain))
            if dim == 2:
                g_ex[0] = 2.+8.*x[0]+5.*x[1]
                g_ex[1] = 3.+5.*x[0]+12.*x[1]
            else:
                g_ex[0] = 2.+6.*x[1]+8.*x[2]+18.*x[0]
                g_ex[1] = 3.+6.*x[0]+7.*x[2]+20.*x[1]
                g_ex[2] = 4.+7.*x[1]+8.*x[0]+22.*x[2]
        return g_ex

    def _getSolution(self, system):
        """returns exact solution"""
        dim = self.domain.getDim()
        x = Solution(self.domain).getX()
        if system:
            u_ex = Vector(0., Solution(self.domain))
            if dim == 2:
                u_ex[0] =  1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
                u_ex[1] = -1.+4.*x[0]+2.*x[1]+1.*x[0]**2+6.*x[1]*x[0]+4.*x[1]**2
            else:
                u_ex[0] = 1.+2.*x[0]+3.*x[1]+4.*x[2]+\
                          6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+\
                          9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
                u_ex[1] = 2.+4.*x[0]+1.*x[1]-6.*x[2]+\
                          3.*x[0]*x[1]+2.*x[1]*x[2]-8.*x[2]*x[0]-\
                          2.*x[0]**2+7.*x[1]**2+5.*x[2]**2
                u_ex[2] = -2.+7.*x[0]+9.*x[1]+2*x[2]-\
                          6.*x[0]*x[1]+8.*x[1]*x[2]+2.*x[2]*x[0]+\
                          2.*x[0]**2+8.*x[1]**2+1.*x[2]**2
        else:
            if dim == 2:
                u_ex = 1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
            else:
                u_ex = 1.+2.*x[0]+3.*x[1]+4.*x[2]+\
                       6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+\
                       9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        return u_ex

    def _setCoefficients(self, pde, system):
        """sets PDE coefficients"""
        super(SimpleSolveOrder2, self)._setCoefficients(pde, system)
        dim = self.domain.getDim()
        if system:
            Y = pde.getCoefficient("Y")
            if dim == 2:
                Y[0] = Y[0]-20.
                Y[1] = Y[1]-10.
            else:
                Y[0] = Y[0]-60.
                Y[1] = Y[1]-20.
                Y[2] = Y[2]-22.
            pde.setValue(Y=Y)
        else:
            if dim == 2:
                pde.setValue(Y=-20.)
            else:
                pde.setValue(Y=-60.)

### BiCGStab + Jacobi

class Test_SimpleSolveFinleyRect_Order1_Paso_BICGSTAB_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_BICGSTAB_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_BICGSTAB_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_BICGSTAB_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### PCG + Jacobi

class Test_SimpleSolveFinleyRect_Order1_Paso_PCG_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_PCG_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_PCG_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_PCG_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### TFQMR + Jacobi

class Test_SimpleSolveFinleyRect_Order1_Paso_TFQMR_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_TFQMR_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_TFQMR_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_TFQMR_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### MINRES + Jacobi

class Test_SimpleSolveFinleyRect_Order1_Paso_MINRES_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_MINRES_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_MINRES_Jacobi(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_MINRES_Jacobi(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### BiCGStab + Gauss-Seidel

class Test_SimpleSolveFinleyRect_Order1_Paso_BICGSTAB_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_BICGSTAB_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_BICGSTAB_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_BICGSTAB_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### PCG + Gauss-Seidel

class Test_SimpleSolveFinleyRect_Order1_Paso_PCG_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_PCG_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_PCG_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_PCG_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### TFQMR + Gauss-Seidel

class Test_SimpleSolveFinleyRect_Order1_Paso_TFQMR_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_TFQMR_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_TFQMR_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_TFQMR_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### MINRES + Gauss-Seidel

class Test_SimpleSolveFinleyRect_Order1_Paso_MINRES_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_MINRES_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_MINRES_GaussSeidel(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_MINRES_GaussSeidel(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### BiCGStab + RILU

class Test_SimpleSolveFinleyRect_Order1_Paso_BICGSTAB_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_BICGSTAB_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_BICGSTAB_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_BICGSTAB_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### PCG + RILU

class Test_SimpleSolveFinleyRect_Order1_Paso_PCG_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_PCG_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_PCG_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_PCG_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### TFQMR + RILU

class Test_SimpleSolveFinleyRect_Order1_Paso_TFQMR_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_TFQMR_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_TFQMR_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_TFQMR_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### MINRES + RILU

class Test_SimpleSolveFinleyRect_Order1_Paso_MINRES_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyRect_Order2_Paso_MINRES_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Rectangle(NE0, NE1, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order1_Paso_MINRES_RILU(SimpleSolveOnPaso):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 1, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

class Test_SimpleSolveFinleyBrick_Order2_Paso_MINRES_RILU(SimpleSolveOrder2):
    def setUp(self):
        self.domain = Brick(NE0, NE1, NE2, 2, optimize=OPTIMIZE)
        self.package = SolverOptions.PASO
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

