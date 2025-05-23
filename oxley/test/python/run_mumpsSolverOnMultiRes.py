
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
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

"""
Test suite for PDE solvers on Oxley multiresolution domains
"""

from test_simplesolve import SimpleSolveTestCase
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

from esys.escript import getMPISizeWorld, hasFeature, sqrt
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import SolverOptions

HAVE_MUMPS = hasFeature('mumps')
skip_amg = True
skip_muelu_long = False #hasFeature("longindex")

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
OPTIMIZE=True
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
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    return m

def test_Rectangle_refine_Point(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refinePoint(x0=0.5,y0=0.5)
    return m

def test_Rectangle_refine_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=0.5)
    return m

def test_Rectangle_refine_Region(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineRegion(x0=0.2,x1=0.2,y0=0.6,y1=0.8)
    return m

def test_Brick_refine_Mesh(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    m.dump("uniform_mesh_ae.silo")
    return m

def test_Brick_refine_Point(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refinePoint(x0=0.55,y0=0.55,z0=0.55)
    m.dump("point_mesh_ae.silo")
    return m

def test_Brick_refine_top_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=DX)
    m.dump("top_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_east_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="right",dx=DX)
    m.dump("east_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_west_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="left",dx=DX)
    m.dump("west_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_bottom_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="bottom",dx=DX)
    m.dump("bottom_boundary_mesh_ae.silo")
    return m

def test_Brick_refine_Region(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m.dump("region_boundary_mesh_ae.silo")
    return m
    return m

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class SimpleSolveOnMumps(SimpleSolveTestCase):
    pass

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Mesh(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Mesh(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Region(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Region(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Region(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Boundary(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Boundary(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Point(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Point(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle_refine_Point(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain                

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Mesh(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Mesh(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Mesh(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Region(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Region(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Region(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Region(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Boundary(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Boundary(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Boundary(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

## direct
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_Direct_Point(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Point(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.DIRECT

    def tearDown(self):
        del self.domain

### BiCGStab + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_Jacobi_Point(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick_refine_Point(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain  

### PCG + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### TFQMR + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_TFQMR_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_TFQMR_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### MINRES + Jacobi
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_MINRES_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_MINRES_Jacobi(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.JACOBI

    def tearDown(self):
        del self.domain

### BiCGStab + Gauss-Seidel
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### PCG + AMG
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

### TFQMR + Gauss-Seidel
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_TFQMR_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_TFQMR_GaussSeidel(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.GAUSS_SEIDEL

    def tearDown(self):
        del self.domain

### MINRES + AMG
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveMultires2D_Mumps_MINRES_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
@unittest.skipIf(skip_muelu_long, "MueLu AMG incompatible with index type long")
@unittest.skipIf(skip_amg, "AMG not available")
class Test_SimpleSolveMultires2D_Mumps_MINRES_AMG(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.AMG

    def tearDown(self):
        del self.domain


### BiCGStab + RILU
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_BICGSTAB_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.BICGSTAB
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### PCG + RILU
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### TFQMR + RILU
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_TFQMR_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_TFQMR_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.TFQMR
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

### MINRES + RILU
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_MINRES_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_MINRES_RILU(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.MINRES
        self.preconditioner = SolverOptions.RILU

    def tearDown(self):
        del self.domainn

### PCG + ILUT
@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_ILUT(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain

@unittest.skipIf(not HAVE_MUMPS, "MUMPS not available")
class Test_SimpleSolveMultires2D_Mumps_PCG_ILUT(SimpleSolveOnMumps):
    def setUp(self):
        self.domain = test_Brick(n0=NE0*NX-1, n1=NE1*NY-1, n2=NE2*NZ-1, d0=NX, d1=NY, d2=NZ)
        self.package = SolverOptions.MUMPS
        self.method = SolverOptions.PCG
        self.preconditioner = SolverOptions.ILUT

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

