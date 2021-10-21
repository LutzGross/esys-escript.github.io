
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

NE=4 # number elements

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


# def Brick(**kwargs):
#     kwargs['n0'] //= 2
#     kwargs['n1'] //= 2
#     kwargs['n2'] //= 2
#     m = MultiResolutionDomain(3, **kwargs)
#     return m.getLevel(1)

class Test_UtilOnOxley_refine_Mesh(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_Mesh(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
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
        self.domain = test_Rectangle_refine_Mesh(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_UtilOnOxley_refine_Point(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_Point(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
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
        self.domain = test_Rectangle_refine_Point(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_UtilOnOxley_refine_Boundary(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_Boundary(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['OXLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_Util_SpatialFunctionsOnOxley2D_refine_Boundary(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = test_Rectangle_refine_Boundary(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_UtilOnOxley_refine_Region(Test_util, Test_symfuncs, Test_util_NaN_funcs):
    def setUp(self):
        self.domain=test_Rectangle_refine_Region(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
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
        self.domain = test_Rectangle_refine_Region(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
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

