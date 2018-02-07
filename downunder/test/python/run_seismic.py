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

import logging
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import numpy as np
import os
import sys
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.escript import *
from esys.weipa import saveSilo

# this is mainly to avoid warning messages
logging.basicConfig(format='%(name)s: %(message)s', level=logging.INFO)

try:
    from esys.ripley import Rectangle as rRect, Brick as rBrick
    HAVE_RIPLEY = True
except ImportError:
    HAVE_RIPLEY = False

try:
    from esys.finley import Rectangle as fRect, Brick as fBrick
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False



try:
    TEST_DATA_ROOT=os.environ['DOWNUNDER_TEST_DATA_ROOT']
except KeyError:
    TEST_DATA_ROOT='ref_data'

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'
    
writeFailMessage = "This feature (SimpleSEGYWriter.write()) depends on obspy,"+\
" which is not installed, see https://github.com/obspy/obspy for install guide"


class TestSeismicTools(unittest.TestCase):
    @unittest.skipIf(getMPISizeWorld() > 1,
            "segywriters don't support multiple ranks")
    def test_segy_writer1(self):
        sw= SimpleSEGYWriter(receiver_group=[0.,1.,2.], source=1., sampling_interval=10*U.msec, text="testing")
        self.assertRaises(ValueError, sw.addRecord, [1])
        for i in range(713):
           # Create some random data.
           data = np.random.ranf(3)
           sw.addRecord(data)
        try:
            sw.write(os.path.join(WORKDIR,"test1.sgy"))
        except Exception as e:
            if str(e) == writeFailMessage:
                raise unittest.SkipTest("obspy not installed")
            raise e

    @unittest.skipIf(getMPISizeWorld() > 1,
            "segywriters don't support multiple ranks")
    def test_segy_writer2(self):
        sw= SimpleSEGYWriter(receiver_group=[(0.,0.),(1.,-1.),(2.,-2)], source=(3,3), sampling_interval=10*U.msec, text="testing")
        self.assertRaises(ValueError, sw.addRecord, [1])
        for i in range(411):
           # Create some random data.
           data = np.random.ranf(3)
           sw.addRecord(data)
        try:
            sw.write(os.path.join(WORKDIR,"test2.sgy"))
        except Exception as e:
            if str(e) == writeFailMessage:
                raise unittest.SkipTest("obspy not installed")
            raise e
        
    def test_ricker(self):
        rw=Ricker(f_dom=40, t_dom=None)

        self.assertLess(abs(rw.getValue(0)), 1e-6)
        self.assertAlmostEqual(rw.getValue(rw.getCenter()), 1. )
        self.assertLess(abs(rw.getAcceleration(0.)), 1.)

    def test_wavebase(self):
        class TestWave(WaveBase):
            def _getAcceleration(self, t, u):
                return -sin(t)
        def check_values(self, tw, t_ref):
            t, u=tw.update(t_ref)
            self.assertLess(abs(t-t_ref), 1e-9)
            self.assertLess(abs(u-sin(t_ref)), 1e-6)
            
        tw=TestWave(dt=0.001, u0=0., v0=1., t0=0.)
        self.assertAlmostEqual(0.001, tw.getTimeStepSize())
        for t in [0.005, 0.007, 0.0071, 0.01, 0.02, 0.5]:
            check_values(self, tw, t)
        self.assertRaises(ValueError, tw.update, t-1)

    def sonicRunner(self, domain, label):
            v_p=1.
            sw = SonicWave(domain, v_p, Ricker(0.5), source_tag='sss')
            u = sw.update(1.)[1]
            self.assertIsInstance(u, Data,
                    "u is not Data instance for %s"%label)
            self.assertEqual(u.getShape(), (),
                    "u is not shape () for %s"%label)
            self.assertEqual(u.getFunctionSpace(), Solution(domain),
                    "functionspace != solution for %s"%label)

    def test_sonicwave2D(self):
        doms = []
        if HAVE_RIPLEY:
            doms.append((rRect, "ripley"))
        if HAVE_FINLEY:
            doms.append((fRect, "finley"))
        for domType, impl in doms:
            domain=domType(5,5, diracPoints=[(0.5,1.)], diracTags=['sss'])
            self.sonicRunner(domain, "%s.Rectangle"%impl)

    def test_sonicwave3D(self):
        doms = []
        if HAVE_RIPLEY:
            doms.append((rBrick, "ripley"))
        if HAVE_FINLEY:
            doms.append((fBrick, "finley"))
        for domType, impl in doms:
            domain=domType(5,5,5, diracPoints=[(0.5,0.5,1.)], diracTags=['sss'])
            self.sonicRunner(domain, "%s.Brick"%impl)
                                  
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
    
