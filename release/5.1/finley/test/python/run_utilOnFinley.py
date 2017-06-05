
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_util import Test_util, Test_Util_SpatialFunctions, \
        Test_Util_SpatialFunctions_noGradOnBoundary, \
        Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_util_NaN_funcs import Test_util_NaN_funcs

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
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

FINLEY_MERGE_ERROR = "merge: more than 1 processor is not supported yet."

NE=4 # number elements, must be even

class Test_UtilOnFinley(Test_util):
   def setUp(self):
       try:
           self.workdir=os.environ['FINLEY_WORKDIR']
       except KeyError:
           self.workdir='.'
       self.domain = Rectangle(NE, NE+1, 2)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_SymFuncsOnFinley(Test_symfuncs,Test_util_NaN_funcs):
   def setUp(self):
       try:
           self.workdir=os.environ['FINLEY_WORKDIR']
       except KeyError:
           self.workdir='.'
       self.domain = Rectangle(NE, NE+1, 2)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_NaNFuncsOnFinley(Test_util_NaN_funcs):
   def setUp(self):
       try:
           self.workdir=os.environ['FINLEY_WORKDIR']
       except KeyError:
           self.workdir='.'
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

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

