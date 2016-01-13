
########################################################
#
# Copyright (c) 2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
import esys.escript
import sys

#Just because we are unit testing these functions does not mean we encourage their use.
#In many cases these tests are merely testing if the call works, not if it does what
#you might expect

class ModuleFnsTestCase(unittest.TestCase):

    def setUp(self):
	pass

    def testGlobalMax(self):
	r=esys.escript.getMPIRankWorld()
	self.assertEqual(esys.escript.getMPISizeWorld()-1,esys.escript.getMPIWorldMax(r))

        
    def testGlobalSum(self):
	r=esys.escript.getMPIRankWorld()
	s=esys.escript.getMPISizeWorld()
	total=s/2.0*(1+s)-s
	self.assertEqual(total,esys.escript.getMPIWorldSum(r))

    def testgetMachinePrecision(self):
	if esys.escript.getMachinePrecision()>1:	#Arbitrary value
		self.fail("Machine precision is not sensible")

    def testMPIBarrier(self):
	esys.escript.MPIBarrierWorld()

    def testgetMaxFloat(self):
	self.assertTrue(esys.escript.getMaxFloat()>1)    #Arbitrary value
	
    def testprintParallelThreadCounts(self):
	esys.escript.printParallelThreadCounts()
	
    def testgetNumberOfThreads(self):
	self.assertTrue(esys.escript.getNumberOfThreads()>=1)
	

    def testgetSvnVersion(self):
	esys.escript.getVersion()
	
if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(ModuleFnsTestCase))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)
