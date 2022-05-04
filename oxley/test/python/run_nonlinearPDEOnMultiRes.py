
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

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

def test_Rectangle_refine_top_Boundary(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=0.5)
    return m

def test_Rectangle_refine_east_Boundary(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="right",dx=0.5)
    return m

def test_Rectangle_refine_west_Boundary(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="left",dx=0.5)
    return m

def test_Rectangle_refine_bottom_Boundary(**kwargs):
    kwargs['n0'] //= 2
    kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="bottom",dx=0.5)
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
# TODO
# @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
# class Test_OxleyNonLinearPDE3D(Test_nlpde):
#    def setUp(self):
#         self.domain = Brick(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
#    def tearDown(self):
#         del self.domain

if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)

