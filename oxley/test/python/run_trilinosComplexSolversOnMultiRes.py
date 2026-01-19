
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


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
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
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_TRILINOS = hasFeature('trilinos')
skip_muelu_long = False #hasFeature("longindex")

# number of elements in the spatial directions
NE0=20
NE1=20
NE2=9

mpiSize=getMPISizeWorld()
# for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#     NX=x
#     NY=mpiSize//x
#     if NX*NY == mpiSize:
#         break

# for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
#     NXb=x[0]
#     NYb=x[1]
#     NZb=mpiSize//(x[0]*x[1])
#     if NXb*NYb*NZb == mpiSize:
#         break

NX=1
NY=1
NZ=1
DX=0.2
ref_level=2

def test_Rectangle_refine_Mesh(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineMesh("uniform")
    m.dump("uniform_mesh_ae.silo")
    return m

def test_Rectangle_refine_Point(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refinePoint(x0=0.55,y0=0.55)
    m.dump("point_mesh_ae.silo")
    return m

def test_Rectangle_refine_top_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="top",dx=DX)
    m.dump("top_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_east_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="right",dx=DX)
    m.dump("east_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_west_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="left",dx=DX)
    m.dump("west_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_bottom_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="bottom",dx=DX)
    m.dump("bottom_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_Region(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m.dump("region_boundary_mesh_ae.silo")
    return m

# TODO
def test_Brick_refine_Mesh(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineMesh("uniform")
    m.dump("uniform_mesh_ae.silo")
    return m

def test_Brick_refine_Point(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refinePoint(x0=0.55,y0=0.55,z0=0.55)
    m.dump("point_mesh_ae.silo")
    return m

def test_Brick_refine_top_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="top",dx=DX)
    m.dump("top_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_east_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="right",dx=DX)
    m.dump("east_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_west_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="left",dx=DX)
    m.dump("west_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_bottom_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineBoundary(boundary="bottom",dx=DX)
    m.dump("bottom_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_Region(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(ref_level)
    m.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m.dump("region_boundary_mesh_ae.silo")
    return m

def Brick(**kwargs):
    m = Brick(**kwargs)
    return m

@unittest.skipIf(not HAVE_TRILINOS, "Trilinos not available")
class ComplexSolveOnTrilinos(ComplexSolveTestCase):
    pass


# TODO
## direct
# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain


# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# TODO
# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain


# class Test_ComplexSolveMultiRes2D_Trilinos_Direct_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.DIRECT

#     def tearDown(self):
#         del self.domain

# ### BiCGStab + Jacobi
# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain
# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

#TODO
# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain
# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# @unittest.skip("fails with Nan during iteration.")
# class Test_ComplexSolveMultiRes2D_Trilinos_BICGSTAB_Jacobi_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.BICGSTAB
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

### GMRES + Jacobi

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# TODO
# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_GMRES_Jacobi_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.GMRES
#         self.preconditioner = SolverOptions.JACOBI
#         self.verbosity = True

#     def tearDown(self):
#         del self.domain

### PCG + Jacobi

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# TODO
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_Jacobi_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.JACOBI

#     def tearDown(self):
#         del self.domain

### PCG + AMG

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# TODO
# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

# @unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_AMG_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.AMG

#     def _setSolverOptions(self, so):
#         so.setTrilinosParameter("number of equations", 2)

#     def tearDown(self):
#         del self.domain

### PCG + ILUT

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1, d0=NX, d1=NY)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# TODO
# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Mesh(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_top_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_east_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_west_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_bottom_Boundary(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Point(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain

# class Test_ComplexSolveMultiRes2D_Trilinos_PCG_ILUT_Region(ComplexSolveOnTrilinos):
#     def setUp(self):
#         self.domain = test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, d0=NX, d1=NY, d2=NZ)
#         self.package = SolverOptions.TRILINOS
#         self.method = SolverOptions.PCG
#         self.preconditioner = SolverOptions.ILUT

#     def tearDown(self):
#         del self.domain


if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)

