
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
Test suite for nonlinearPDEs on Oxley
"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_nonLinearPDE import Test_nlpde
from esys.escript import getMPISizeWorld
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick

mpiSize = getMPISizeWorld()

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


# def Brick(**kwargs):
#     m = MultiResolutionDomain(3, **kwargs)
#     return m.getLevel(1)

class Test_OxleyNonLinearPDE2D_Mesh(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle_Mesh(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_Point(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle_Point(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_top_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle_top_Boundary(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_east_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle_east_Boundary(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_west_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle_west_Boundary(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_bottom_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle_bottom_Boundary(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain


class Test_OxleyNonLinearPDE2D_Region(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle_Region(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain                

class Test_OxleyNonLinearPDE2D_Mesh(Test_nlpde):
   def setUp(self):
        self.domain = test_Brick_refine_Mesh(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_Point(Test_nlpde):
   def setUp(self):
        self.domain = test_Brick_refine_Point(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_top_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Brick_refine_top_Boundary(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_east_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Brick_refine_east_Boundary(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_west_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Brick_refine_west_Boundary(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

class Test_OxleyNonLinearPDE2D_bottom_Boundary(Test_nlpde):
   def setUp(self):
        self.domain = test_Brick_refine_bottom_Boundary(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain


class Test_OxleyNonLinearPDE2D_Region(Test_nlpde):
   def setUp(self):
        self.domain = test_Brick_refine_Region(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)

