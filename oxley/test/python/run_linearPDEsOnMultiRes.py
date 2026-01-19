
########################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test suite for the linearPDE and pdetools on Oxley

:remark:

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
# from test_linearPDEs import Test_Poisson, Test_LinearPDE, Test_TransportPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping
from test_linearPDEs import Test_Poisson, Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_3Do1
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import *
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick

try:
     OXLEY_TEST_DATA=os.environ['OXLEY_TEST_DATA']
except KeyError:
     OXLEY_TEST_DATA='.'

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
DX=0.2

def test_Rectangle_refine_Mesh(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    m.dump("uniform_mesh_ae.silo")
    return m

def test_Rectangle_refine_Point(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refinePoint(x0=0.55,y0=0.55)
    m.dump("point_mesh_ae.silo")
    return m

def test_Rectangle_refine_top_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=DX)
    m.dump("top_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_east_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="right",dx=DX)
    m.dump("east_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_west_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="left",dx=DX)
    m.dump("west_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_bottom_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="bottom",dx=DX)
    m.dump("bottom_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_Region(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m.dump("region_boundary_mesh_ae.silo")
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

# class Test_LinearPDEOnOxleyRect(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1, Test_TransportPDE):
@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRect_Mesh(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRect_Point(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_Point(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRect_top_Boundary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, l0=10., l1=10., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRect_east_Boundary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, l0=10., l1=10., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRect_west_Boundary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, l0=10., l1=10., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRect_bottom_Boundary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, l0=10., l1=10., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRect_Region(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_Region(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

# TODO
# class Test_LinearPDEOnOxleyBrick_Mesh(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1
#     def tearDown(self):
#         del self.domain

# class Test_LinearPDEOnOxleyRBrickPoint(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1
#     def tearDown(self):
#         del self.domain

# class Test_LinearPDEOnOxleyRect_topBrickndary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=10., l1=10., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1
#     def tearDown(self):
#         del self.domain

# class Test_LinearPDEOnOxleyRect_eastBrickndary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=10., l1=10., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1
#     def tearDown(self):
#         del self.domain

# class Test_LinearPDEOnOxleyRect_westBrickndary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=10., l1=10., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1
#     def tearDown(self):
#         del self.domain

# class Test_LinearPDEOnOxleyRect_bottomBrickndary(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=10., l1=10., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1
#     def tearDown(self):
#         del self.domain

# class Test_LinearPDEOnOxleyReBrickegion(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1
#     def tearDown(self):
#         del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley_Mesh(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley_Point(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_Point(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley_top_Boundary(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley_east_Boundary(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley_west_Boundary(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley_bottom_Boundary(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley_Region(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=test_Rectangle_refine_Region(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain

# TODO
# class Test_PoissonOnOxleyBrick_Mesh(Test_Poisson):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#     def tearDown(self):
#         del self.domain

# class Test_PoissonOnOxleyBrick_Point(Test_Poisson):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#     def tearDown(self):
#         del self.domain

# class Test_PoissonOnOxleyBrick_top_Boundary(Test_Poisson):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_top_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#     def tearDown(self):
#         del self.domain

# class Test_PoissonOnOxleyBrick_east_Boundary(Test_Poisson):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_east_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#     def tearDown(self):
#         del self.domain

# class Test_PoissonOnOxleyBrick_west_Boundary(Test_Poisson):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_west_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#     def tearDown(self):
#         del self.domain

# class Test_PoissonOnOxleyBrick_bottom_Boundary(Test_Poisson):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#     def tearDown(self):
#         del self.domain

# class Test_PoissonOnOxleyBrick_Region(Test_Poisson):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#             NX=x
#             NY=mpiSize//x
#             NZ=x
#             if NX*NY == mpiSize:
#                 break
#         self.domain=test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#     def tearDown(self):
#         del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

