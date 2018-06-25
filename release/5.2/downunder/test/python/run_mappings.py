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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2012-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys
from esys.escript import *
from esys.downunder.mappings import *
import numpy

try:
    from esys.ripley import Rectangle
    HAVE_RIPLEY = True
except ImportError:
    HAVE_RIPLEY = False

class Test_LinearMapping(unittest.TestCase):
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

@unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
class Test_AcousticVelocityMapping(unittest.TestCase):
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
        STEP=0.5e-4
        
        DS=mp.getDerivative((m0, m1))
        s0=mp.getValue((m0, m1))
        s1=mp.getValue((m0+STEP, m1))
        DS_ref = (s1-s0)/STEP
        self.assertLess( Lsup(DS[0,0]-DS_ref[0]), 1e-2 * Lsup(DS_ref[0]))
        self.assertLess( Lsup(DS[1,0]-DS_ref[1]), 1e-2 * Lsup(DS_ref[1]))

        s1=mp.getValue((m0, m1+STEP))
        DS_ref= (s1-s0)/STEP
        self.assertTrue( Lsup(DS[1,1]-DS_ref[1]), 1e-2 * Lsup(DS_ref[1]))
        self.assertLess( Lsup(DS[0,1]-DS_ref[0]), 1e-2 * Lsup(DS_ref[0]))
        # get inverse
        s0_ref=1/(V0*complex(1,-1/(2*Q0)))**2
        m0=mp.getInverse( (s0_ref.real,s0_ref.imag))
        self.assertLess( Lsup(m0[0]), 1e-14)
        self.assertLess( Lsup(m0[1]), 1e-14)

class Test_BoundedRangeMapping(unittest.TestCase):
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

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

