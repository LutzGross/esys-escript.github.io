
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
from test_util import Test_util, Test_Util_SpatialFunctions, \
        Test_Util_SpatialFunctions_noGradOnBoundary, \
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_types import Test_addition_types
from test_util_interpolation import Test_Util_Interpolation_Dirac
from test_util_integrals import Test_Util_Integration_Dirac
from test_util_NaN_funcs import Test_util_NaN_funcs
# from test_copyWithData import Test_copyWithMask

from esys.escript import FunctionOnBoundary, getMPISizeWorld, HAVE_SYMBOLS
from esys.finley import Rectangle, Brick, JoinFaces, ReadMesh
import os

if HAVE_SYMBOLS:
    from test_symfuncs import Test_symfuncs
else:
    @unittest.skip("Skipping symbolic tests since sympy is not available")
    class Test_symfuncs:
        pass

try:
     ESCRIPT_WORKDIR=os.environ['ESCRIPT_WORKDIR']
except KeyError:
     ESCRIPT_WORKDIR='.'

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

FINLEY_MERGE_ERROR = "merge: more than 1 processor is not supported yet."

NE=4 # number elements, must be even

#TODO: What is this testing? Is this still needed?
class Test_Data_Addition(Test_addition_types):
    def setUp(self):
        self.domain = Rectangle(NE, NE+1, 2)
        self.workdir=ESCRIPT_WORKDIR

        self.functionspace = FunctionOnBoundary(self.domain)
    def tearDown(self):
       del self.domain


class Test_UtilOnFinley(Test_util,Test_util_NaN_funcs):
   def setUp(self):
       self.workdir = ESCRIPT_WORKDIR
       self.domain = Rectangle(NE, NE + 1, 2)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain

# class Test_SymFuncsOnFinley(Test_symfuncs):
#    def setUp(self):
#        try:
#            self.workdir=os.environ['FINLEY_WORKDIR']
#        except KeyError:
#            self.workdir='.'
#        self.domain = Rectangle(NE, NE+1, 2)
#        self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
#    def tearDown(self):
#        del self.functionspace
#        del self.domain

class Test_NaNFuncsOnFinley(Test_util_NaN_funcs):
   def setUp(self):
       self.workdir = ESCRIPT_WORKDIR
       self.domain = Rectangle(NE, NE+1, 2)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet2DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet2DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=2
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_order2.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet2DMacro(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_2D_macro.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet3DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet3DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=2
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_order2.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet3DMacro(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet_3D_macro.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Rectangle(n0=NE,n1=NE,order=1,useElementsOnFace=0)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=2
        self.domain = Rectangle(n0=NE,n1=NE,order=2,useElementsOnFace=0)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DMacro(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Rectangle(n0=NE,n1=NE,order=-1,useElementsOnFace=0)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Brick(n0=NE,n1=NE,n2=NE,order=1,useElementsOnFace=0)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=2
        self.domain = Brick(n0=NE,n1=NE,n2=NE,order=2,useElementsOnFace=0)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DMacro(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Brick(n0=NE,n1=NE,n2=NE,order=-1,useElementsOnFace=0)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder1withContact(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=1
        d1 = Rectangle(n0=NE//2+1,n1=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2 = Rectangle(n0=NE//2,n1=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder2withContact(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=2
        d1 = Rectangle(n0=NE//2,n1=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2 = Rectangle(n0=NE//2,n1=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder1withContact(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=1
        d1 = Brick(n0=NE//2+1,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2 = Brick(n0=NE//2,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder2withContact(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=2
        d1 = Brick(n0=NE//2+1,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2 = Brick(n0=NE//2,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder1useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Rectangle(n0=NE//2+1,n1=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Rectangle(n0=NE//2,n1=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder2useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Rectangle(n0=NE//2+1,n1=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Rectangle(n0=NE//2,n1=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder1useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Brick(n0=NE//2,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Brick(n0=NE//2+1,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder2useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Brick(n0=NE//2,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Brick(n0=NE//2+1,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        if getMPISizeWorld() > 1:
            with self.assertRaises(NotImplementedError) as pkg:
                self.domain = JoinFaces([d1,d2],optimize=False)
            e = pkg.exception
            if FINLEY_MERGE_ERROR not in str(e):
                raise e
            raise unittest.SkipTest(FINLEY_MERGE_ERROR)
        else:
            self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_2D_Point_Order1_Interpolation(Test_Util_Interpolation_Dirac):
    def setUp(self):
        Stations = [ (0.,0.), (1.,0), (0,1), (1,1) ]
        StationsTags = ["A1", "A2", "A3", "A4" ]
        self.positions=Stations
        self.taglist=StationsTags
        self.domain=Rectangle(n0=5,n1=5, diracPoints=Stations, diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.positions, self.taglist

class Test_3D_Point_Order1_Interpolation(Test_Util_Interpolation_Dirac):
    def setUp(self):
        Stations = [ (0.,0.,0.), (1.,0,0.), (0,1,0.), (1,1,0.), (0.,0.,1.), (1.,0,1.), (0,1,1.), (1,1,1.) ]
        StationsTags = ["A1", "A2", "A3", "A4","A5", "A6", "A7", "A8"  ]
        self.positions=Stations
        self.taglist=StationsTags
        self.domain=Brick(n0=5,n1=5,n2=5,diracPoints=Stations,diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.positions, self.taglist

class Test_2D_Point_Order2_Interpolation(Test_Util_Interpolation_Dirac):
    def setUp(self):
        Stations = [ (0.,0.), (1.,0), (0,1), (1,1) ]
        StationsTags = ["A1", "A2", "A3", "A4" ]
        self.positions=Stations
        self.taglist=StationsTags
        self.domain=Rectangle(n0=5,n1=5, order=2, diracPoints=Stations, diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.positions, self.taglist

class Test_3D_Point_Order2_Interpolation(Test_Util_Interpolation_Dirac):
    def setUp(self):
        Stations = [ (0.,0.,0.), (1.,0,0.), (0,1,0.), (1,1,0.), (0.,0.,1.), (1.,0,1.), (0,1,1.), (1,1,1.) ]
        StationsTags = ["A1", "A2", "A3", "A4","A5", "A6", "A7", "A8"  ]
        self.positions=Stations
        self.taglist=StationsTags
        self.domain=Brick(n0=5,n1=5,n2=5,order=2, diracPoints=Stations,diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.positions, self.taglist

class Test_2D_Point_Order1_Integration(Test_Util_Integration_Dirac):
    def setUp(self):
        Stations = [ (0.,0.), (1.,0), (0,1), (1,1) ]
        StationsTags = ["A1", "A2", "A3", "A4" ]
        self.taglist=StationsTags
        self.domain=Rectangle(n0=5,n1=5, diracPoints=Stations, diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist

class Test_3D_Point_Order1_Integration(Test_Util_Integration_Dirac):
    def setUp(self):
        Stations = [ (0.,0.,0.), (1.,0,0.), (0,1,0.), (1,1,0.), (0.,0.,1.), (1.,0,1.), (0,1,1.), (1,1,1.) ]
        StationsTags = ["A1", "A2", "A3", "A4","A5", "A6", "A7", "A8"  ]
        self.taglist=StationsTags
        self.domain=Brick(n0=5,n1=5,n2=5,diracPoints=Stations,diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist
class Test_2D_Point_Order2_Integration(Test_Util_Integration_Dirac):
    def setUp(self):
        Stations = [ (0.,0.), (1.,0), (0,1), (1,1) ]
        StationsTags = ["A1", "A2", "A3", "A4" ]
        self.taglist=StationsTags
        self.domain=Rectangle(n0=5,n1=5, order=2, diracPoints=Stations, diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist

class Test_3D_Point_Order2_Integration(Test_Util_Integration_Dirac):
    def setUp(self):
        Stations = [ (0.,0.,0.), (1.,0,0.), (0,1,0.), (1,1,0.), (0.,0.,1.), (1.,0,1.), (0,1,1.), (1,1,1.) ]
        StationsTags = ["A1", "A2", "A3", "A4","A5", "A6", "A7", "A8"  ]
        self.taglist=StationsTags
        self.domain=Brick(n0=5,n1=5,n2=5,order=2, diracPoints=Stations,diracTags=StationsTags)
    def tearDown(self):
        del self.domain, self.taglist

# TODO
# class Test_copyWithMask_rectangle(Test_copyWithMask):
#     def setUp(self):
#         self.domain=Rectangle(n0=5,n1=7)
#     def tearDown(self):
#         del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
