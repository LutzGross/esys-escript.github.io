
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import logging
import unittest
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
    TEST_DATA_ROOT=os.environ['DOWNUNDER_TEST_DATA_ROOT']
except KeyError:
    TEST_DATA_ROOT='ref_data'

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'
    



class TestSeismicTools(unittest.TestCase):
    def test_segy_writer1(self):
        sw= SimpleSEGYWriter(receiver_group=[0.,1.,2.], source=1., sampling_interval=10*U.msec, text="testing")
        self.assertRaises(ValueError, sw.addRecord, [1])
        for i in range(713):
           # Create some random data.
           data = np.random.ranf(3)
           sw.addRecord(data)
        try:
            sw.write(os.path.join(WORKDIR,"test1.sgy"))
        except Exception:
            raise unittest.SkipTest("obspy may not be installed")
        
    def test_segy_writer2(self):
        sw= SimpleSEGYWriter(receiver_group=[(0.,0.),(1.,-1.),(2.,-2)], source=(3,3), sampling_interval=10*U.msec, text="testing")
        self.assertRaises(ValueError, sw.addRecord, [1])
        for i in range(411):
           # Create some random data.
           data = np.random.ranf(3)
           sw.addRecord(data)
        try:
            sw.write(os.path.join(WORKDIR,"test2.sgy"))
        except Exception:
            raise unittest.SkipTest("obspy may not be installed")
        
    def test_ricker(self):
        rw=Ricker(f_dom=40, t_dom=None)

        self.assertTrue(abs(rw.getValue(0)) < 1e-6)
        self.assertAlmostEquals(rw.getValue(rw.getCenter()), 1. )
        self.assertTrue(abs(rw.getAcceleration(0.)) < 1.)

    def test_wavebase(self):
        class TestWave(WaveBase):
            def _getAcceleration(self, t, u):
                return -sin(t)
        tw=TestWave(dt=0.001, u0=0., v0=1., t0=0.)
        self.assertAlmostEquals(0.001, tw.getTimeStepSize())

        t_ref=0.005
        t, u=tw.update(t_ref)
        self.assertTrue(abs(t-t_ref) < 1e-9)
        self.assertTrue(abs(u-sin(t_ref)) < 1e-6)

        t_ref=0.007
        t, u=tw.update(t_ref)
        self.assertTrue(abs(t-t_ref) < 1e-9)
        self.assertTrue(abs(u-sin(t_ref)) < 1e-6)

        t_ref=0.0071
        t, u=tw.update(t_ref)
        self.assertTrue(abs(t-t_ref) < 1e-9)
        self.assertTrue(abs(u-sin(t_ref)) < 1e-6)

        t_ref=0.01
        t, u=tw.update(t_ref)
        self.assertTrue(abs(t-t_ref) < 1e-9)
        self.assertTrue(abs(u-sin(t_ref)) < 1e-6)

        t_ref=0.02
        t, u=tw.update(t_ref)
        self.assertTrue(abs(t-t_ref) < 1e-9)
        self.assertTrue(abs(u-sin(t_ref)) < 1e-6)

        t_ref=0.5
        t, u=tw.update(t_ref)
        self.assertTrue(abs(t-t_ref) < 1e-9)
        self.assertTrue(abs(u-sin(t_ref)) < 1e-6)

        self.assertRaises(ValueError, tw.update, t_ref-1)

    def test_sonicwave2D(self):
        from esys.finley import Rectangle

        domain=Rectangle(5,5, diracPoints=[(0.5,1.)], diracTags=['sss'])
        v_p=1.

        sw=SonicWave( domain, v_p, Ricker(0.5), source_tag='sss')
        u=sw.update(1.)[1]
        self.assertTrue(isinstance(u,Data))
        self.assertEqual(u.getShape(), ())
        self.assertEqual(u.getFunctionSpace(), Solution(domain))

        from esys.ripley import Rectangle

        domain=Rectangle(5,5, diracPoints=[(0.5,1.)], diracTags=['sss'])
        v_p=1.

        sw=SonicWave( domain, v_p, Ricker(0.5), source_tag='sss')
        u=sw.update(1.)[1]
        self.assertTrue(isinstance(u,Data))
        self.assertEqual(u.getShape(), ())
        self.assertEqual(u.getFunctionSpace(), Solution(domain))
         
    def test_sonicwave3D(self):
        from esys.finley import Brick

        domain=Brick(5,5,5, diracPoints=[(0.5,0.5,1.)], diracTags=['sss'])
        v_p=1.

        sw=SonicWave( domain, v_p, Ricker(0.5), source_tag='sss')
        u=sw.update(1.)[1]
        self.assertTrue(isinstance(u,Data))
        self.assertEqual(u.getShape(), ())
        self.assertEqual(u.getFunctionSpace(), Solution(domain))
        
        from esys.ripley import Brick

        domain=Brick(5,5,5, diracPoints=[(0.5,0.5,1.)], diracTags=['sss'])
        v_p=1.

        sw=SonicWave( domain, v_p, Ricker(0.5), source_tag='sss')
        u=sw.update(1.)[1]
        self.assertTrue(isinstance(u,Data))
        self.assertEqual(u.getShape(), ())
        self.assertEqual(u.getFunctionSpace(), Solution(domain))
                 
                                  
if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSeismicTools))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)
    
