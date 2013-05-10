
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
from esys.downunder.mappings import *

class TestLinearMapping(unittest.TestCase):
    def setUp(self):
        self.alpha = 4.2
        self.p0 = 1.75
        self.mapping=LinearMapping(self.alpha, self.p0)

    def test_invalid_alpha(self):
        self.assertRaises(Exception, LinearMapping, 0)

    def test_value(self):
        v=self.mapping.getValue(1.23)
        ref=self.alpha*1.23+self.p0
        self.assertAlmostEqual(v, ref)

    def test_derivative(self):
        v=self.mapping.getDerivative(1.23)
        ref=self.alpha
        self.assertAlmostEqual(v, ref)

    def test_inverse(self):
        v=self.mapping.getInverse(1.23)
        ref=(1.23-self.p0)/self.alpha
        self.assertAlmostEqual(v, ref)

class TestBoundedRangeMapping(unittest.TestCase):
    def setUp(self):
        self.alpha=4.2
        self.mapping=BoundedRangeMapping(1.5, 3.5)

    def test_invalid_range(self):
        self.assertRaises(ValueError, BoundedRangeMapping, 1, 0)

    def test_min_value(self):
        v=self.mapping.getValue(-10000)
        ref=1.5
        self.assertAlmostEqual(v, ref)

    def test_max_value(self):
        v=self.mapping.getValue(10000)
        ref=3.5
        self.assertAlmostEqual(v, ref)

    def test_zero_value(self):
        v=self.mapping.getValue(0)
        ref=2.5
        self.assertAlmostEqual(v, ref)

    def test_min_derivative(self):
        v=self.mapping.getDerivative(-10000)
        ref=0.
        self.assertAlmostEqual(v, ref)

    def test_max_derivative(self):
        v=self.mapping.getDerivative(10000)
        ref=0.
        self.assertAlmostEqual(v, ref)

    def test_zero_derivative(self):
        v=self.mapping.getDerivative(0)
        ref=1.
        self.assertAlmostEqual(v, ref)

    def test_inverse_toolow(self):
        self.assertRaises(ValueError, self.mapping.getInverse, 1.5)

    def test_inverse_toohigh(self):
        self.assertRaises(ValueError, self.mapping.getInverse, 3.5)

    def test_inverse(self):
        v=self.mapping.getInverse(2.5)
        ref=0.
        self.assertAlmostEqual(v, ref)

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestLinearMapping))
    suite.addTest(unittest.makeSuite(TestBoundedRangeMapping))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

