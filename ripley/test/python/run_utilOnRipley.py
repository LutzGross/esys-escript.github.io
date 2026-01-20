
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
from test_util_interpolation import Test_Util_Interpolation_Dirac
from test_util_integrals import Test_Util_Integration_Dirac
# from test_util_interpolation import Test_Util_Point_Data_Interpolation
from test_symfuncs import Test_symfuncs
from esys.escript import *
from esys.ripley import Rectangle, Brick

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

class Test_UtilOnRipley(Test_util,Test_util_NaN_funcs):
    def setUp(self):
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['RIPLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain

class Test_SymFuncsOnRipley(Test_symfuncs):
    def setUp(self):
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
        try:
            self.workdir=os.environ['RIPLEY_WORKDIR']
        except KeyError:
            self.workdir='.'

    def tearDown(self):
        del self.functionspace
        del self.domain
        
class Test_Util_SpatialFunctionsOnRipley2D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnRipley3D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Brick(n0=NE*NXb-1, n1=NE*NYb-1, n2=NE*NZb-1, l0=1., l1=1., l2=1., d0=NXb, d1=NYb, d2=NZb)
    def tearDown(self):
        del self.order
        del self.domain

class Test_2D_Point_Data_Interpolation(Test_Util_Interpolation_Dirac):
    def setUp(self):
        Stations = [ (0.,0.), (1.,0), (0,1), (1,1) ]
        StationsTags = ["A1", "A2", "A3", "A4" ]
        self.taglist=StationsTags
        self.positions=Stations
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY, diracPoints=Stations, diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist, self.positions

class Test_3D_Point_Data_Interpolation(Test_Util_Interpolation_Dirac):
    def setUp(self):
        Stations = [ (0.,0.,0.), (1.,0,0.), (0,1,0.), (1,1,0.), (0.,0.,1.), (1.,0,1.), (0,1,1.), (1,1,1.) ]
        StationsTags = ["A1", "A2", "A3", "A4","A5", "A6", "A7", "A8"  ]
        self.taglist=StationsTags
        self.positions=Stations
        self.domain=Brick(n0=NE*NXb-1, n1=NE*NYb-1, n2=NE*NZb-1, l0=1., l1=1., l2=1., d0=NXb, d1=NYb, d2=NZb ,diracPoints=Stations,diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist, self.positions

class Test_2D_Point_Data_Integration(Test_Util_Integration_Dirac):
    def setUp(self):
        Stations = [ (0.,0.), (1.,0), (0,1), (1,1) ]
        StationsTags = ["A1", "A2", "A3", "A4" ]
        self.taglist=StationsTags
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY, diracPoints=Stations, diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist

class Test_3D_Point_Data_Integration(Test_Util_Integration_Dirac):
    def setUp(self):
        Stations = [ (0.,0.,0.), (1.,0,0.), (0,1,0.), (1,1,0.), (0.,0.,1.), (1.,0,1.), (0,1,1.), (1,1,1.) ]
        StationsTags = ["A1", "A2", "A3", "A4","A5", "A6", "A7", "A8"  ]
        self.taglist=StationsTags
        self.domain=Brick(n0=NE*NXb-1, n1=NE*NYb-1, n2=NE*NZb-1, l0=1., l1=1., l2=1., d0=NXb, d1=NYb, d2=NZb ,diracPoints=Stations,diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist
        



if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
