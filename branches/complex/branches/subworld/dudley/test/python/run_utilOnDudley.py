
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

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_util import Test_util as Test_util
from test_util import Test_Util_SpatialFunctions, Test_Util_SpatialFunctions_noGradOnBoundary_noContact

from esys.escript import *
from esys.dudley import Rectangle,Brick,ReadMesh
import sys
import os

if HAVE_SYMBOLS:
    from test_symfuncs import Test_symfuncs
else:
    print("Skipping symbolic tests since sympy is not available")
    class Test_symfuncs:
        pass

try:
     DUDLEY_TEST_DATA=os.environ['DUDLEY_TEST_DATA']
except KeyError:
     DUDLEY_TEST_DATA='.'

DUDLEY_TEST_MESH_PATH=os.path.join(DUDLEY_TEST_DATA,"data_meshes")


NE=4 # number elements, must be even

class Test_UtilOnDudley(Test_util,Test_symfuncs):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,1)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_Util_SpatialFunctionsOnDudleyTet2DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain


class Test_Util_SpatialFunctionsOnDudleyTet3DOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnDudleyRectOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Rectangle(n0=NE,n1=NE,order=1)
    def tearDown(self):
        del self.order
        del self.domain


class Test_Util_SpatialFunctionsOnDudleyBrickOrder1(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Brick(n0=NE,n1=NE,n2=NE,order=1)
    def tearDown(self):
        del self.order
        del self.domain



if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

