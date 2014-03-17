
##############################################################################
#
# Copyright (c) 2012-2014 by University of Queensland
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

__copyright__="""Copyright (c) 2012-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
import sys
from esys.escript import *
from esys.downunder.mappings import *
from esys.ripley import Rectangle
import numpy

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
        
class TestAcousticVelocityMapping(unittest.TestCase):
    def setUp(self):
        self.fs=Function(Rectangle(40,40))
    def tearDown(self):
        del self.fs

    def test_case1(self):
        Q0=100.
        V0=3000
        QQ=Scalar(Q0, self.fs)
        VV=Scalar(V0, self.fs)
        mp=AcousticVelocityMapping(VV, QQ)
        
        # get Value test
        s0=mp.getValue( (0,0))
        s0_ref=1/(V0*complex(1,-1/(2*Q0)))**2
        self.assertLess( Lsup(s0[0]-s0_ref.real), 1e-8 * abs(s0_ref.real))
        self.assertLess( Lsup(s0[1]-s0_ref.imag), 1e-8 * abs(s0_ref.imag))
        
        M0=numpy.log(s0_ref)
        s1=mp.getValue( (1.,0))
        s1_ref=numpy.exp(complex(1.,0)+M0)
        self.assertLess( Lsup(s1[0]-s1_ref.real), 1e-8 * abs(s1_ref.real))
        self.assertLess( Lsup(s1[1]-s1_ref.imag), 1e-8 * abs(s1_ref.imag))
        
        s2=mp.getValue( (0,1.))
        s2_ref=numpy.exp(complex(0,1,)+M0)
        self.assertLess( Lsup(s2[0]-s2_ref.real), 1e-8 * abs(s2_ref.real))
        self.assertLess( Lsup(s2[1]-s2_ref.imag), 1e-8 * abs(s2_ref.imag))
        
        # get derivative
        m0=0.5
        m1=0.01
        STEP=1e-3
        
        DS=mp.getDerivative((m0, m1))
        s0=mp.getValue((m0, m1))
        s1=mp.getValue((m0+STEP, m1))
        DS_ref = (s1-s0)/STEP
        self.assertLess( Lsup(DS[0,0]-DS_ref[0]), 1e-8 * Lsup(DS_ref[0]))
        self.assertLess( Lsup(DS[1,0]-DS_ref[1]), 1e-8 * Lsup(DS_ref[1]))

        s1=mp.getValue((m0, m1+STEP))
        DS_ref= (s1-s0)/STEP
        self.assertLess( Lsup(DS[0,1]-DS_ref[0]), 1e-8 * Lsup(DS_ref[0]))
        self.assertTrue( Lsup(DS[1,1]-DS_ref[1]), 1e-8 * Lsup(DS_ref[1]))
        # get inverse
        s0_ref=1/(V0*complex(1,-1/(2*Q0)))**2
        m0=mp.getInverse( (s0_ref.real,s0_ref.imag))
        self.assertLess( Lsup(m0[0]), 1e-14)
        self.assertLess( Lsup(m0[1]), 1e-14)
        
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
    suite.addTest(unittest.makeSuite(TestAcousticVelocityMapping))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

