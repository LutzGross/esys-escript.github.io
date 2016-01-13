
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
from test_util import Test_util as Test_util
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary, Test_Util_SpatialFunctions_noGradOnBoundary_noContact

from esys.escript import *
from esys.finley import Rectangle,Brick,JoinFaces,ReadMesh
import sys
import os

if HAVE_SYMBOLS:
    from test_symfuncs import Test_symfuncs
else:
    print("Skipping symbolic tests since sympy is not available")
    class Test_symfuncs:
        pass

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")


NE=4 # number elements, must be even

class Test_UtilOnFinley(Test_util,Test_symfuncs):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
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
        d1 = Rectangle(n0=NE/2+1,n1=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder2withContact(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=2
        d1 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder1withContact(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=1
        d1 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder2withContact(Test_Util_SpatialFunctions_noGradOnBoundary):
    def setUp(self):
        self.order=2
        d1 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=0)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder1useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Rectangle(n0=NE/2+1,n1=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex2DOrder2useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Rectangle(n0=NE/2+1,n1=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder1useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=1
        d1 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnFinleyHex3DOrder2useElementsOnFacewithContact(Test_Util_SpatialFunctions):
    def setUp(self):
        self.order=2
        d1 = Brick(n0=NE/2,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2 = Brick(n0=NE/2+1,n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
        d2.setX(d2.getX()+[0.5,0.,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   if True:
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DOrder1))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DOrder2))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DMacro))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex3DOrder1))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex3DOrder2))
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex3DMacro))
   else:
      suite.addTest(unittest.makeSuite(Test_Util_SpatialFunctionsOnFinleyHex2DOrder1useElementsOnFacewithContact))
      pass
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

