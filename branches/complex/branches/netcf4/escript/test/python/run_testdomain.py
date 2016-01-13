
##############################################################################
#
# Copyright (c) 2012-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2012-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
import sys
from esys.escript import *

class TestDomainTests(unittest.TestCase):
    TOL=EPSILON*100.
    def testreduction(self):
        dom = getTestDomainFunctionSpace(4,2,2).getDomain()
        dx=dom.getX()
        msk_arg=whereNegative(dx[0]-0.5)
        arg=msk_arg*numpy.array([-0.5, -0.84])+(1.-msk_arg)*numpy.array([-0.72, 0.9])
        res=maxval(arg)
        msk_ref=whereNegative(dx[0]-0.5)
        ref=msk_ref*(-0.5)+(1.-msk_ref)*0.9
        self.assertTrue(Lsup(res-ref) <= self.TOL, "ReductionOnTestDomain Failed")

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDomainTests))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)
