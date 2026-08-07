
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

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_util import Test_util
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_symfuncs import Test_symfuncs
from esys.escript import *
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick
from esys.oxley import RefinementQueue2D, RefinementQueue3D

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
    _f = RefinementQueue2D()
    _f.setRefinementLevel(1)
    _f.refineUniform()
    m = _f.apply(m)
    m.dump("uniform_mesh_ae.silo")
    return m

def test_Rectangle_refine_Point(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    _f = RefinementQueue2D()
    _f.setRefinementLevel(1)
    _f.refinePoint(x0=0.55,y0=0.55)
    m = _f.apply(m)
    m.dump("point_mesh_ae.silo")
    return m

def test_Rectangle_refine_top_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    _f = RefinementQueue2D()
    _f.setRefinementLevel(1)
    _f.refineBorder(border="top",dx=DX)
    m = _f.apply(m)
    m.dump("top_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_east_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    _f = RefinementQueue2D()
    _f.setRefinementLevel(1)
    _f.refineBorder(border="right",dx=DX)
    m = _f.apply(m)
    m.dump("east_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_west_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    _f = RefinementQueue2D()
    _f.setRefinementLevel(1)
    _f.refineBorder(border="left",dx=DX)
    m = _f.apply(m)
    m.dump("west_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_bottom_Boundary(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    _f = RefinementQueue2D()
    _f.setRefinementLevel(1)
    _f.refineBorder(border="bottom",dx=DX)
    m = _f.apply(m)
    m.dump("bottom_boundary_mesh_ae.silo")
    return m

def test_Rectangle_refine_Region(**kwargs):
    # kwargs['n0'] //= 2
    # kwargs['n1'] //= 2
    m = Rectangle(**kwargs)
    _f = RefinementQueue2D()
    _f.setRefinementLevel(1)
    _f.refineRegion(x0=0.2,x1=0.6,y0=0.6,y1=0.8)
    m = _f.apply(m)
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

