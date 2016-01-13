
##############################################################################
#
# Copyright (c) 2009-2016 by The University of Queensland
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2009-2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
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
        if esys.escript.getMachinePrecision()>1:        #Arbitrary value
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
        
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
