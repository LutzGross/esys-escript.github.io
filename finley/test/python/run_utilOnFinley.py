# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import unittest
from test_util import Test_util as Test_util
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary, Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_symbols import Test_symbols

from esys.escript import FunctionOnBoundary, saveVTK
from esys.finley import Rectangle,Brick,JoinFaces,ReadMesh
import sys
import os

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=FINLEY_TEST_DATA+"/data_meshes/"


NE=4 # number elements, must be even

class Test_UtilOnFinley(Test_util,Test_symbols):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet2DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_2D_order1.fly")
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet2DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=2
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_2D_order2.fly")
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet3DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_3D_order1.fly")
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyTet3DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=2
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_3D_order2.fly")
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=1
        d1 = Rectangle(n0=NE/2+1,n1=NE,l0=0.5,order=1)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=1)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=2
        d1 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=2)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=2)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=1
        d1 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=1)
        d2 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=1)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder2(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=2
        d1 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=2)
        d2 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=2)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder1useElementsOnFace(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Rectangle(n0=NE/2+1,n1=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder2useElementsOnFace(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Rectangle(n0=NE/2+1,n1=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder1useElementsOnFace(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder2useElementsOnFace(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2])
    def tearDown(self):
        del self.order
        del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   # suite.addTest(Test_UtilOnFinleyHex("test_generalTensorTransposedProduct_Symbol_rank1_taggedData_rank3_offset1"))
   suite.addTest(unittest.makeSuite(Test_UtilOnFinley))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DOrder1))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DOrder2))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex3DOrder1))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex3DOrder2))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyTet2DOrder1))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyTet2DOrder2))
   # suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyTet3DOrder1))
   # suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyTet3DOrder2))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DOrder1useElementsOnFace))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DOrder2useElementsOnFace))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex3DOrder1useElementsOnFace))
   suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex3DOrder2useElementsOnFace))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
   
