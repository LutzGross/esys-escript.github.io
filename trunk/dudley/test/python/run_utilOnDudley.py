
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
from test_util import Test_Util_SpatialFunctions_noGradOnBoundary_noContact

from esys.escript import FunctionOnBoundary, HAVE_SYMBOLS
from esys.dudley import Rectangle,Brick,ReadMesh
import os

if HAVE_SYMBOLS:
    from test_symfuncs import Test_symfuncs
else:
    @unittest.skip("Skipping symbolic tests since sympy is not available")
    class Test_symfuncs:
        pass

try:
     DUDLEY_TEST_DATA=os.environ['DUDLEY_TEST_DATA']
except KeyError:
     DUDLEY_TEST_DATA='.'

try:
     DUDLEY_WORKDIR=os.environ['DUDLEY_WORKDIR']
except KeyError:
     DUDLEY_WORKDIR='.'

DUDLEY_TEST_MESH_PATH=os.path.join(DUDLEY_TEST_DATA,"data_meshes")


NE=4 # number elements, must be even

class Test_UtilOnDudley(Test_util):
   def setUp(self):
       self.domain = Rectangle(NE,NE+1)
       self.functionspace = FunctionOnBoundary(self.domain) # due to a bug in escript python needs to hold a reference to the domain
       self.workdir=DUDLEY_WORKDIR

   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_SymFuncsOnDudley(Test_symfuncs):
   def setUp(self):
       self.domain = Rectangle(NE,NE+1)
       self.functionspace = FunctionOnBoundary(self.domain)
       self.workdir=DUDLEY_WORKDIR

   def tearDown(self):
       del self.functionspace
       del self.domain

class Test_Util_SpatialFunctionsOnDudleyTet2D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_2D_order1.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain


class Test_Util_SpatialFunctionsOnDudleyTet3D(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = ReadMesh(os.path.join(DUDLEY_TEST_MESH_PATH,"tet_3D_order1.fly"),optimize=False)
    def tearDown(self):
        del self.order
        del self.domain

class Test_Util_SpatialFunctionsOnDudleyRect(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Rectangle(n0=NE,n1=NE)
    def tearDown(self):
        del self.order
        del self.domain


class Test_Util_SpatialFunctionsOnDudleyBrick(Test_Util_SpatialFunctions_noGradOnBoundary_noContact):
    def setUp(self):
        self.order=1
        self.domain = Brick(n0=NE,n1=NE,n2=NE)
    def tearDown(self):
        del self.order
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

