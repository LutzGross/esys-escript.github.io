#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
Test suite for input and output of meshes and data objects

@remark:

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Ken Steube, k.steube@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision: 859 $"
__date__="$Date: 2006-09-26 12:19:18 +1000 (Tue, 26 Sep 2006) $"

import unittest, sys

from esys.escript import *
from esys.finley import Rectangle, Brick, LoadMesh, ReadMesh, ReadMeshMPI

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=FINLEY_TEST_DATA+"/data_meshes/"

SOLVER_TOL=1.e-8
REL_TOL=1.e-6

class InputOutput(unittest.TestCase):

     def domainsEqual(self, m1, m2, nft=20):
        self.failUnless(m1.getDim() == m2.getDim(), "Dimensions differ")
        self.failUnless(m1.getNumDataPointsGlobal() == m2.getNumDataPointsGlobal(), "Num data points differ")
        for tagName in m1.showTagNames().split(", "):
          self.failUnless(m2.isValidTagName(tagName), "m1 has a tag '%s' not present in m2" % tagName)
        for tagName in m2.showTagNames().split(", "):
          self.failUnless(m1.isValidTagName(tagName), "m2 has a tag '%s' not present in m1" % tagName)
          self.failUnless(m1.getTag(tagName) == m2.getTag(tagName), "values of tag '%s' differ" % tagName)
        for fs in ["Solution", "ReducedSolution", "Function", "ReducedFunction", "FunctionOnBoundary", "ReducedFunctionOnBoundary", "ContinuousFunction", "ReducedContinuousFunction"]:
          fs1 = eval("%s(m1)" % fs)
          fs2 = eval("%s(m2)" % fs)
          x1 = fs1.getX()
          x2 = fs2.getX()
          for n in range(1, nft+1):
            integ1 = integrate(sin(n*x1))
            integ2 = integrate(sin(n*x2))
            self.failUnless(Lsup(abs(integ1-integ2)) <= REL_TOL, "integrals for n=%d differ" % n)
        return True

     # Does optimize=True change Rectangle for order=1?
     def test_Rectangle_optimize_order1(self):
	mydomain1 = Rectangle(n0=7, n1=11, order=1, l0=1., l1=1., optimize=False)
	mydomain2 = Rectangle(n0=7, n1=11, order=1, l0=1., l1=1., optimize=True)
        self.domainsEqual(mydomain1, mydomain2)

     # Does optimize=True change Rectangle for order=2?
     def test_Rectangle_optimize_order2(self):
	mydomain1 = Rectangle(n0=7, n1=11, order=2, l0=1., l1=1., optimize=False)
	mydomain2 = Rectangle(n0=7, n1=11, order=2, l0=1., l1=1., optimize=True)
        self.domainsEqual(mydomain1, mydomain2)

     # Does optimize=True change Brick for order=1?
     def test_Brick_optimize_order1(self):
	mydomain1 = Brick(n0=7, n1=11, n2=5, order=1, l0=1., l1=1., l2=1., optimize=False)
	mydomain2 = Brick(n0=7, n1=11, n2=5, order=1, l0=1., l1=1., l2=1., optimize=True)
        self.domainsEqual(mydomain1, mydomain2)

     # Does optimize=True change Brick for order=2?
     def test_Brick_optimize_order2(self):
	mydomain1 = Brick(n0=7, n1=11, n2=5, order=2, l0=1., l1=1., l2=1., optimize=False)
	mydomain2 = Brick(n0=7, n1=11, n2=5, order=2, l0=1., l1=1., l2=1., optimize=True)
        self.domainsEqual(mydomain1, mydomain2)

     def test_mesh_dump_to_NetCDF_rectangle(self):
	if loadIsConfigured():
	  mydomain1 = Rectangle(n0=7, n1=11, order=1, l0=1., l1=1., optimize=False)
	  mydomain1.dump("tt.mesh.nc")
	  mydomain2=LoadMesh("tt.mesh.nc")
          self.domainsEqual(mydomain1, mydomain2)

     def DISABLED_test_mesh_dump_to_NetCDF_brick(self):
	if loadIsConfigured():
	  mydomain1 = Brick(n0=7, n1=11, n2=5, order=2, l0=1., l1=1., l2=1., optimize=False)
	  mydomain1.dump("tt.mesh.nc")
	  mydomain2=LoadMesh("tt.mesh.nc")
          self.domainsEqual(mydomain1, mydomain2)

     def test_mesh_read_rectangle_from_finley_file(self):
	mydomain1 = Rectangle(n0=8, n1=10, order=1, l0=1., l1=1., optimize=False)
        mydomain2 = ReadMeshMPI(FINLEY_TEST_MESH_PATH+"rectangle_8x10.fly")
        self.domainsEqual(mydomain1, mydomain2)

     def test_mesh_read_brick_from_finley_file(self):
        mydomain1 = Brick(n0=8, n1=10, n2=12, order=1, l0=1., l1=1., l2=1., optimize=False)
        mydomain2 = ReadMeshMPI(FINLEY_TEST_MESH_PATH+"brick_8x10x12.fly")
        self.domainsEqual(mydomain1, mydomain2)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(InputOutput))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

