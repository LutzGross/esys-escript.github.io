
##############################################################################
#
# Copyright (c) 2012-2018 by The University of Queensland
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


__copyright__="""Copyright (c) 2012-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

from esys.escript import *
import numpy
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys


class TestDomainTests(unittest.TestCase):
    """
    escript has a (relatively) trivial domain implementatiobn (called TestDomain) 
    to be used for testing escriptcore code without needing a full domain 
    implementation such as finley.
    
    This set of tests checks that domain.
    """
    EPSILON=1e-8
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

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
