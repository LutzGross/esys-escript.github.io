
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
Test suite for PDE solvers on Oxley
"""

from test_simplesolve import ComplexSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_TRILINOS = hasFeature('trilinos')
skip_muelu_long = False #hasFeature("longindex")

# number of elements in the spatial directions
NE0=10
NE1=10
NE2=9
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

def test_Rectangle_refine_Mesh(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    return m

def test_Rectangle_refine_Point(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refinePoint(x0=0.5,y0=0.5)
    return m

def test_Rectangle_refine_Boundary(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=0.5)
    return m

def test_Rectangle_refine_Region(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineRegion(x0=0.2,x1=0.2,y0=0.6,y1=0.8)
    return m

# def Brick(**kwargs):
#     m = MultiResolutionDomain(3, **kwargs)
#     return m.getLevel(1)

@unittest.skipIf(not HAVE_TRILINOS, "Trilinos not available")
class ComplexSolveOnTrilinos(ComplexSolveTestCase):
    pass


## direct
class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Mesh(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Point(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Boundary(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Region(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

# TODO
# class Test_ComplexSolveMultiRes3D_Trilinos_Direct(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

### BiCGStab + Jacobi
@unittest.skip("fails with Nan during iteration.")
class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Mesh(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

@unittest.skip("fails with Nan during iteration.")
class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Boundary(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

@unittest.skip("fails with Nan during iteration.")
class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Point(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

@unittest.skip("fails with Nan during iteration.")
class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Region(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

# TODO
# class Test_ComplexSolveMultiRes3D_Trilinos_BICGSTAB_Jacobi(ComplexSolveOnTrilinos):
#     SOLVER_TOL = 1.e-9
#     def setUp(self):
#         self.domain = Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

### GMRES + Jacobi

class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Mesh(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.GMRES
        self.preconditioner = SolverOptions.JACOBI
        self.verbosity = True

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Point(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.GMRES
        self.preconditioner = SolverOptions.JACOBI
        self.verbosity = True

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Boundary(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.GMRES
        self.preconditioner = SolverOptions.JACOBI
        self.verbosity = True

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Region(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.GMRES
        self.preconditioner = SolverOptions.JACOBI
        self.verbosity = True

    def tearDown(self):
        del self.domain
# TODO
# class Test_ComplexSolveMultiRes3D_Trilinos_GMRES_Jacobi(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

### PCG + Jacobi

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Mesh(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Point(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Boundary(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Region(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

# TODO
# class Test_ComplexSolveMultiRes3D_Trilinos_PCG_Jacobi(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

### PCG + AMG

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Mesh(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def _setSolverOptions(self, so):
        so.setTrilinosParameter("number of equations", 2)

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Boundary(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def _setSolverOptions(self, so):
        so.setTrilinosParameter("number of equations", 2)

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Point(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def _setSolverOptions(self, so):
        so.setTrilinosParameter("number of equations", 2)

    def tearDown(self):
        del self.domain

@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Region(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def _setSolverOptions(self, so):
        so.setTrilinosParameter("number of equations", 2)

    def tearDown(self):
        del self.domain

# TODO
# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes3D_Trilinos_PCG_AMG(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 3)

#     def tearDown(self):
#         del self.domain

### PCG + ILUT

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Mesh(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Boundary(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Point(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Region(ComplexSolveOnTrilinos):
    def setUp(self):
        self.domain = test_Rectangle_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.TRILINOS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

# TODO
# class Test_ComplexSolveMultiRes3D_Trilinos_PCG_ILUT(ComplexSolveOnTrilinos):
#     SOLVER_TOL = 1.e-9
#     def setUp(self):
#         self.domain = Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain


if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)

