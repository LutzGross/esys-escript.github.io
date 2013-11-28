
##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
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
        for i in xrange(713):
           # Create some random data.
           data = np.random.ranf(3)
           sw.addRecord(data)
        sw.write(os.path.join(WORKDIR,"test1.sgy"))
        
     def test_segy_writer2(self):
        sw= SimpleSEGYWriter(receiver_group=[(0.,0.),(1.,-1.),(2.,-2)], source=(3,3), sampling_interval=10*U.msec, text="testing")
        self.assertRaises(ValueError, sw.addRecord, [1])
        for i in xrange(411):
           # Create some random data.
           data = np.random.ranf(3)
           sw.addRecord(data)
        sw.write(os.path.join(WORKDIR,"test2.sgy"))



if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSeismicTools))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)
    
