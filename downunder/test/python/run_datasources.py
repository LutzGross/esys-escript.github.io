
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
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
from esys.escript import inf,sup,saveDataCSV
from esys.downunder.datasources import *

# this is mainly to avoid warning messages
logger=logging.getLogger('inv')
logger.setLevel(logging.INFO)
handler=logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)

try:
    TEST_DATA_ROOT=os.environ['DOWNUNDER_TEST_DATA_ROOT']
except KeyError:
    TEST_DATA_ROOT='.'

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'


ERS_DATA = os.path.join(TEST_DATA_ROOT, 'ermapper_test.ers')
ERS_REF = os.path.join(TEST_DATA_ROOT, 'ermapper_test.csv')
ERS_NULL = -99999 * 1e-6
ERS_SIZE = [20,15]
ERS_ORIGIN = [309097.0, 6319002.0]
NC_DATA = os.path.join(TEST_DATA_ROOT, 'netcdf_test.nc')
NC_REF = os.path.join(TEST_DATA_ROOT, 'netcdf_test.csv')
NC_NULL = 0.
NC_SIZE = [20,15]
NC_ORIGIN = [309097.0, 6319002.0]
VMIN=-10000.
VMAX=10000
NE_V=15
ALT=0.
PAD_X=7
PAD_Y=9

class TestERSDataSource(unittest.TestCase):
    def test_ers_with_padding(self):
        source = ERSDataSource(headerfile=ERS_DATA, vertical_extents=(VMIN,VMAX,NE_V), alt_of_data=ALT)
        source.setPadding(PAD_X,PAD_Y)
        dom=source.getDomain()
        g,s=source.getGravityAndStdDev()

        outfn=os.path.join(WORKDIR, '_ersdata.csv')
        saveDataCSV(outfn, g=g[2], s=s)

        X0,NP,DX=source.getDataExtents()
        V0,NV,DV=source.getVerticalExtents()

        # check metadata
        self.assertEqual(NP, ERS_SIZE, msg="Wrong number of data points")
        # this test only works if gdal is available
        try:
            import osgeo.osr
            for i in xrange(len(ERS_ORIGIN)):
                self.assertAlmostEqual(X0[i], ERS_ORIGIN[i], msg="Data origin wrong")
        except ImportError:
            print("Skipping test of data origin since gdal is not installed.")

        # check data
        nx=NP[0]+2*PAD_X
        ny=NP[1]+2*PAD_Y
        nz=NE_V
        z_data=int(np.round((ALT-V0)/DV)-1)

        ref=np.genfromtxt(ERS_REF, delimiter=',', dtype=float)
        g_ref=ref[:,0].reshape((NP[1],NP[0]))
        s_ref=ref[:,1].reshape((NP[1],NP[0]))

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=float)
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)/(nx*ny)
        g_out=out[:,0].reshape(nz,ny,nx)
        s_out=out[:,1].reshape(nz,ny,nx)
        self.assertAlmostEqual(np.abs(
            g_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]-g_ref).max(),
            0., msg="Difference in gravity data area")

        self.assertAlmostEqual(np.abs(
            s_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]-s_ref).max(),
            0., msg="Difference in error data area")

        # overwrite data -> should only be padding value left
        g_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]=ERS_NULL
        self.assertAlmostEqual(np.abs(g_out-ERS_NULL).max(), 0.,
                msg="Wrong values in padding area")

class TestNetCDFDataSource(unittest.TestCase):
    def test_cdf_with_padding(self):
        source = NetCDFDataSource(gravfile=NC_DATA, vertical_extents=(VMIN,VMAX,NE_V), alt_of_data=ALT)
        source.setPadding(PAD_X,PAD_Y)
        dom=source.getDomain()
        g,s=source.getGravityAndStdDev()

        outfn=os.path.join(WORKDIR, '_ncdata.csv')
        saveDataCSV(outfn, g=g[2], s=s)

        X0,NP,DX=source.getDataExtents()
        V0,NV,DV=source.getVerticalExtents()

        # check metadata
        self.assertEqual(NP, NC_SIZE, msg="Wrong number of data points")
        # this only works if gdal is available
        #self.assertAlmostEqual(X0, NC_ORIGIN, msg="Data origin wrong")

        # check data
        nx=NP[0]+2*PAD_X
        ny=NP[1]+2*PAD_Y
        nz=NE_V
        z_data=int(np.round((ALT-V0)/DV)-1)

        ref=np.genfromtxt(NC_REF, delimiter=',', dtype=float)
        g_ref=ref[:,0].reshape((NP[1],NP[0]))
        s_ref=ref[:,1].reshape((NP[1],NP[0]))

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=float)
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)/(nx*ny)
        g_out=out[:,0].reshape(nz,ny,nx)
        s_out=out[:,1].reshape(nz,ny,nx)

        self.assertAlmostEqual(np.abs(
            g_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]-g_ref).max(),
            0., msg="Difference in gravity data area")

        self.assertAlmostEqual(np.abs(
            s_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]-s_ref).max(),
            0., msg="Difference in error data area")

        # overwrite data -> should only be padding value left
        g_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]=NC_NULL
        self.assertAlmostEqual(np.abs(g_out-NC_NULL).max(), 0.,
                msg="Wrong values in padding area")

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestERSDataSource))
    if 'NetCDFDataSource' in dir():
        suite.addTest(unittest.makeSuite(TestNetCDFDataSource))
    else:
        print("Skipping netCDF data source test since netCDF is not installed")
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

