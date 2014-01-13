
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
import numpy as np
import os
import sys
import unittest
from esys.escript import inf,sup,saveDataCSV,getMPISizeWorld
from esys.downunder.datasources import *
from esys.downunder.domainbuilder import DomainBuilder
from esys.downunder import WGS84ReferenceSystem


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


NC_DATA1 = os.path.join(TEST_DATA_ROOT, 'zone51.nc')
NC_DATA2 = os.path.join(TEST_DATA_ROOT, 'zone52.nc')

@unittest.skipIf('NetCdfData' not in dir(), 'netCDF not available')
class TestDomainBuilderWithNetCdf(unittest.TestCase):

    # append custom messages instead of overwriting originals
    longMessage = True

    def test_add_garbage(self):
        db=DomainBuilder()
        self.assertRaises(TypeError, db.addSource, 42)

    def test_mixing_utm_zones(self):
        source1 = NetCdfData(DataSource.GRAVITY, NC_DATA1, scale_factor=1.)
        source2 = NetCdfData(DataSource.GRAVITY, NC_DATA2, scale_factor=1.)
        domainbuilder=DomainBuilder()
        domainbuilder.addSource(source1)
        self.assertRaises(ValueError, domainbuilder.addSource, source2)

    def test_add_source_after_domain_built(self):
        db=DomainBuilder()
        source1a = NetCdfData(DataSource.GRAVITY, NC_DATA1, scale_factor=1.)
        db.addSource(source1a)
        _=db.getDomain()
        source1b = NetCdfData(DataSource.GRAVITY, NC_DATA1, scale_factor=2.)
        self.assertRaises(Exception, db.addSource, source1b)

    def test_cartesian_domain(self):
        db=DomainBuilder()
        source1a = NetCdfData(DataSource.GRAVITY, NC_DATA1, scale_factor=1.)
        db.addSource(source1a)
        db.setVerticalExtents(depth=20000., air_layer=30000., num_cells=10)
        dom=db.getDomain()

    def test_geodetic_domain(self):
        COORDINATES=WGS84ReferenceSystem()
        db=DomainBuilder(reference_system=COORDINATES)
        source1a = NetCdfData(DataSource.GRAVITY, NC_DATA1, scale_factor=1., reference_system=COORDINATES)
        db.addSource(source1a)
        db.setVerticalExtents(depth=20000., air_layer=30000., num_cells=10)
        dom=db.getDomain()
        x=dom.getX()
        self.assertAlmostEqual(inf(x[0]), 120.2, delta=0.001, msg="phi range wrong")
        self.assertAlmostEqual(inf(x[1]), -29.2 , delta=0.0001, msg="lambda range wrong")
        self.assertAlmostEqual(inf(x[2]), -0.2, msg="h range wrong"+str(x[2]))
        # Cannot check upper bounds of coordinates with more than 1 rank
        # because ripley may adjust internally.
        if getMPISizeWorld()==1:
            self.assertAlmostEqual(sup(x[0]), 120.3, delta=0.001, msg="phi range wrong")
            self.assertAlmostEqual(sup(x[1]), -29.1, delta=0.0001, msg="lambda range wrong")
            self.assertAlmostEqual(sup(x[2]), 0.3, msg="h range wrong: "+str(x[2]))
        

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDomainBuilderWithNetCdf))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

