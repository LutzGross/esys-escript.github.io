
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
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import LameEquation

# from run_customAssemblersOnOxley import OxleyLameAssemblerTestBase, OxleyWaveAssemblerTestBase, Ricker
from run_customAssemblersOnOxley import OxleyLameAssemblerTestBase #, Ricker

NE0=20
NE1=20
NE2=10

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
DY=0.2

def test_Rectangle_refine_Mesh(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    m.dump("uniform_mesh.silo")
    return m

def test_Rectangle_refine_Point(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(2)
    m.refinePoint(x0=DX,y0=DY)
    m.dump("point_mesh.silo")
    return m

def test_Rectangle_refine_top_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=DX)
    m.dump("top_boundary_mesh.silo")
    return m

def test_Rectangle_refine_east_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="right",dx=DX)
    m.dump("east_boundary_mesh.silo")
    return m

def test_Rectangle_refine_west_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="left",dx=DX)
    m.dump("west_boundary_mesh.silo")
    return m

def test_Rectangle_refine_bottom_Boundary(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="bottom",dx=DX)
    m.dump("bottom_boundary_mesh.silo")
    return m

def test_Rectangle_refine_Region(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m.dump("region_boundary_mesh.silo")
    return m

def test_Brick_refine_Mesh(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    m.dump("uniform_mesh.silo")
    return m

def test_Brick_refine_Point(**kwargs):
    # m = Brick(**kwargs)
    m = Brick(n0=10,n1=20,n2=20,l0=1.0,l1=2.0,l2=2.0)
    m.setRefinementLevel(1)
    m.refinePoint(x0=0.55,y0=0.55,z0=0.55)
    m.dump("point_mesh.silo")
    return m

def test_Brick_refine_top_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=DX)
    m.dump("top_boundary_mesh.silo")
    return m

def test_Brick_refine_east_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="east",dx=DX)
    m.dump("east_boundary_mesh.silo")
    return m

def test_Brick_refine_west_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="west",dx=DX)
    m.dump("west_boundary_mesh.silo")
    return m

def test_Brick_refine_bottom_Boundary(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="bottom",dx=DX)
    m.dump("bottom_boundary_mesh.silo")
    return m

def test_Brick_refine_Region(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m.dump("region_boundary_mesh.silo")
    return m

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors - see issue #118")
class Test_OxleyLameAssemblers2D_Mesh(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1)

    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors - see issue #118")
class Test_OxleyLameAssemblers2D_Point(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1)

    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors - see issue #118")
class Test_OxleyLameAssemblers2D_top_Boundary(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1,l0=10,l1=10)

    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors - see issue #118")
class Test_OxleyLameAssemblers2D_east_Boundary(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1,l0=10,l1=10)

    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors - see issue #118")
class Test_OxleyLameAssemblers2D_west_Boundary(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1,l0=10,l1=10)

    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors - see issue #118")
class Test_OxleyLameAssemblers2D_bottom_Boundary(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1,l0=10,l1=10)

    def tearDown(self):
        del self.domain

@unittest.skip("Oxley multi-resolution meshes have SystemMatrixPattern errors - see issue #118")
class Test_OxleyLameAssemblers2D_Region(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1)

    def tearDown(self):
        del self.domain

# TODO
# class Test_OxleyLameAssemblers3D_Mesh(OxleyLameAssemblerTestBase):
#     def setUp(self):
#         self.domain = test_Brick_refine_Mesh(n0=NE0, n1=NE1, n2=NE2)

#     def tearDown(self):
#         del self.domain

# class Test_OxleyLameAssemblers3D_Point(OxleyLameAssemblerTestBase):
#     def setUp(self):
#         self.domain = test_Brick_refine_Point(n0=NE0, n1=NE1, n2=NE2)

#     def tearDown(self):
#         del self.domain

# class Test_OxleyLameAssemblers3D_top_Boundary(OxleyLameAssemblerTestBase):
#     def setUp(self):
#         self.domain = test_Brick_refine_top_Boundary(n0=NE0, n1=NE1,l0=10,l1=10, n2=NE2)

#     def tearDown(self):
#         del self.domain

# class Test_OxleyLameAssemblers3D_east_Boundary(OxleyLameAssemblerTestBase):
#     def setUp(self):
#         self.domain = test_Brick_refine_east_Boundary(n0=NE0, n1=NE1,l0=10,l1=10, n2=NE2)

#     def tearDown(self):
#         del self.domain

# class Test_OxleyLameAssemblers3D_west_Boundary(OxleyLameAssemblerTestBase):
#     def setUp(self):
#         self.domain = test_Brick_refine_west_Boundary(n0=NE0, n1=NE1,l0=10,l1=10, n2=NE2)

#     def tearDown(self):
#         del self.domain

# class Test_OxleyLameAssemblers3D_bottom_Boundary(OxleyLameAssemblerTestBase):
#     def setUp(self):
#         self.domain = test_Brick_refine_bottom_Boundary(n0=NE0, n1=NE1,l0=10,l1=10, n2=NE2)

#     def tearDown(self):
#         del self.domain

# class Test_OxleyLameAssemblers3D_Region(OxleyLameAssemblerTestBase):
#     def setUp(self):
#         self.domain = test_Brick_refine_Region(n0=NE0, n1=NE1, n2=NE2)

#     def tearDown(self):
#         del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

