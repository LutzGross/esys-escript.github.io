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

def domainsEqual(m1, m2, nft=10):
   if m1.getDim() != m2.getDim():
     print "domainsEqual: different dimensions %d, %d" % (m1.getDim(), m2.getDim())
     return False
   if m1.getNumDataPointsGlobal() != m2.getNumDataPointsGlobal():
     print "domainsEqual: different num data points: %d, %d" % (m1.getNumDataPointsGlobal(), m2.getNumDataPointsGlobal())
     return False
   for tagName in m1.showTagNames().split(", "):
     if not m2.isValidTagName(tagName):
       print "domainsEqual: m1 has a tag '%s' not present in m2" % tagName
       return False
   for tagName in m2.showTagNames().split(", "):
     if not m1.isValidTagName(tagName):
       print "domainsEqual: m2 has a tag '%s' not present in m1" % tagName
       return False
     if m1.getTag(tagName) != m2.getTag(tagName):
       print "domainsEqual: values of tag '%s' differ" % tagName
       return False
   for fs in ["Solution", "ReducedSolution", "Function", "ReducedFunction", "FunctionOnBoundary", "ReducedFunctionOnBoundary", "ContinuousFunction", "ReducedContinuousFunction"]:
     fs1 = eval("%s(m1)" % fs)
     fs2 = eval("%s(m2)" % fs)
     x1 = fs1.getX()
     x2 = fs2.getX()
     for n in range(1, nft+1):
       integ1 = integrate(sin(n*x1))
       integ2 = integrate(sin(n*x2))
       if Lsup(abs(integ1-integ2)) > REL_TOL:
	 print "domainsEqual: integrals for n=%d differ" % n
         return False
   return True

class InputOutput(unittest.TestCase):

     # Does optimize=True change Rectangle for order=1?
     def test_Rectangle_optimize_order1(self):
	mydomain1 = Rectangle(n0=7, n1=11, order=1, l0=1., l1=1., optimize=False)
	mydomain2 = Rectangle(n0=7, n1=11, order=1, l0=1., l1=1., optimize=True)
        self.failUnless(domainsEqual(mydomain1, mydomain2), "Domains differ")

     # Does optimize=True change Rectangle for order=2?
     def test_Rectangle_optimize_order2(self):
	mydomain1 = Rectangle(n0=7, n1=11, order=2, l0=1., l1=1., optimize=False)
	mydomain2 = Rectangle(n0=7, n1=11, order=2, l0=1., l1=1., optimize=True)
        self.failUnless(domainsEqual(mydomain1, mydomain2), "Domains differ")

     # Does optimize=True change Brick for order=1?
     def test_Brick_optimize_order1(self):
	mydomain1 = Brick(n0=7, n1=11, n2=5, order=1, l0=1., l1=1., l2=1., optimize=False)
	mydomain2 = Brick(n0=7, n1=11, n2=5, order=1, l0=1., l1=1., l2=1., optimize=True)
        self.failUnless(domainsEqual(mydomain1, mydomain2), "Domains differ")

     # Does optimize=True change Brick for order=2?
     def test_Brick_optimize_order2(self):
	mydomain1 = Brick(n0=7, n1=11, n2=5, order=2, l0=1., l1=1., l2=1., optimize=False)
	mydomain2 = Brick(n0=7, n1=11, n2=5, order=2, l0=1., l1=1., l2=1., optimize=True)
        self.failUnless(domainsEqual(mydomain1, mydomain2), "Domains differ")

     def DISABLED_test_NetCDF(self):
	if loadIsConfigured():
	  mydomain1 = Rectangle(n0=2, n1=3, order=1, l0=1., l1=1.,
	    integrationOrder=-1, reducedIntegrationOrder=-1, periodic0=False, periodic1=False,
	    useElementsOnFace=False, useFullElementOrder=False, optimize=False)
	  mydomain1.dump("tt.mesh.nc")
	  mydomain2=LoadMesh("tt.mesh.nc")
          self.failUnless(domainsEqual(mydomain1, mydomain2), "Domains differ")

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(InputOutput))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

