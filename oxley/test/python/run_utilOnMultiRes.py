
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
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_util import Test_util
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_symfuncs import Test_symfuncs
from esys.escript import *
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick

if HAVE_SYMBOLS:
    from test_symfuncs import Test_symfuncs
else:
    print("Skipping symbolic tests since sympy is not available")
    class Test_symfuncs:
        pass

from test_util_NaN_funcs import Test_util_NaN_funcs

NE0=20 # number elements
NE1=20 # number elements

mpiSize=getMPISizeWorld()
NX=1
NY=1
DX=0.2

def test_Rectangle_refine_Mesh(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    m.dump("uniform_mesh_ae.silo")
    return m

def test_Rectangle_refine_Point(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refinePoint(x0=0.55,y0=0.55)
    m.dump("point_mesh_ae.silo")
    return m

def test_Rectangle_refine_top_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="top",dx=DX)
    m.dump("top_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_east_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="right",dx=DX)
    m.dump("east_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_west_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="left",dx=DX)
    m.dump("west_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_bottom_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineBoundary(boundary="bottom",dx=DX)
    m.dump("bottom_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_Region(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m.dump("region_boundary_mesh_ae.silo")
    return m


# def Brick(**kwargs):
#     kwargs['n0'] //= 2
#     kwargs['n1'] //= 2
#     kwargs['n2'] //= 2
#     m = MultiResolutionDomain(3, **kwargs)
#     return m.getLevel(1)

class Test_UtilOnOxley_refine_Mesh(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_Mesh(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_Mesh(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_UtilOnOxley_refine_Point(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_Point(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_Point(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_Point(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_UtilOnOxley_refine_top_Boundary(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_UtilOnOxley_refine_east_Boundary(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_UtilOnOxley_refine_west_Boundary(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_UtilOnOxley_refine_bottom_Boundary(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_top_Boundary(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_top_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_east_Boundary(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_east_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_west_Boundary(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_west_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_bottom_Boundary(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_bottom_Boundary(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain


class Test_UtilOnOxley_refine_Region(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_Region(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_Region(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_Region(n0=NE0, n1=NE1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

# TODO
# @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
# class Test_Util_SpatialFunctionsOnOxley3D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
#     def setUp(self):
#         self.order=1
#         self.domain = Brick(n0=NE*NXb-1, n1=NE*NYb-1, n2=NE*NZb-1, l0=1., l1=1., l2=1., d0=NXb, d1=NYb, d2=NZb)
#     def tearDown(self):
#         del self.order
#         del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

