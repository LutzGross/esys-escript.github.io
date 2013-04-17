
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

try:
    import pyproj
    haveProj=True
except ImportError:
    haveProj=False
mpisize = getMPISizeWorld()

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


ERS_DATA = os.path.join(TEST_DATA_ROOT, 'ermapper_test.ers')
ERS_REF = os.path.join(TEST_DATA_ROOT, 'ermapper_test.csv')
ERS_NULL = -99999 * 1e-6
ERS_SIZE = [20,15]
ERS_ORIGIN = [309241.0, 6318655.0]
NC_DATA = os.path.join(TEST_DATA_ROOT, 'netcdf_test.nc')
NC_REF = os.path.join(TEST_DATA_ROOT, 'netcdf_test.csv')
NC_NULL = 0.
NC_SIZE = [20,15]
NC_ORIGIN = [403320.91466610413, 6414860.942530109]
NUMPY_NULL = -123.4
VMIN=-10000.
VMAX=10000.
NE_V=15
ALT=0.
PAD_X=3
PAD_Y=2

class TestNumpyData(unittest.TestCase):
    def test_numpy_argument_check(self):
        # invalid data type
        self.assertRaises(ValueError, NumpyData, '_mydatatype_', [1,2])
        # invalid shape of data
        self.assertRaises(ValueError, NumpyData, DataSource.GRAVITY, 42)
        # invalid shape of data
        self.assertRaises(ValueError, NumpyData, DataSource.GRAVITY, np.zeros((2,2,2,2)))
        # invalid shape of error
        self.assertRaises(ValueError, NumpyData, DataSource.GRAVITY, [1,2], [1,2,3])
        # invalid shape of length
        self.assertRaises(ValueError, NumpyData, DataSource.GRAVITY, [1,2], [1,2], [2,3,2])

    @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
    def test_numpy_data_1d(self):
        DIM=1
        testdata = np.arange(20)
        error = 1.*np.ones(testdata.shape)
        source = NumpyData(DataSource.GRAVITY, testdata, null_value=NUMPY_NULL)
        X0,NP,DX=source.getDataExtents()
        for i in range(DIM):
            self.assertAlmostEqual(X0[i], 0., msg="Data origin wrong")
            self.assertEqual(NP[i], testdata.shape[DIM-i-1], msg="Wrong number of data points")
            self.assertAlmostEqual(DX[i], 1000./testdata.shape[DIM-i-1], msg="Wrong cell size")

        domainbuilder=DomainBuilder(dim=2)
        domainbuilder.addSource(source)
        domainbuilder.setVerticalExtents(depth=-VMIN, air_layer=VMAX, num_cells=NE_V)
        domainbuilder.setElementPadding(PAD_X)
        dom=domainbuilder.getDomain()
        g,s=domainbuilder.getGravitySurveys()[0]

        outfn=os.path.join(WORKDIR, '_npdata1d.csv')
        saveDataCSV(outfn, g=g, s=s)

        DV=(VMAX-VMIN)/NE_V

        # check data
        nx=NP[0]+2*PAD_X
        nz=NE_V
        z_data=int(np.round((ALT-VMIN)/DV)-1)

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=float)
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)/nx
        g_out=out[:,0].reshape(nz,nx)
        s_out=out[:,1].reshape(nz,nx)
        self.assertAlmostEqual(np.abs(
            g_out[z_data, PAD_X:PAD_X+NP[0]]-testdata).max(),
            0., msg="Difference in gravity data area")

        self.assertAlmostEqual(np.abs(
            s_out[z_data, PAD_X:PAD_X+NP[0]]-error).max(),
            0., msg="Difference in error data area")

        # overwrite data -> should only be padding value left
        g_out[z_data, PAD_X:PAD_X+NP[0]]=NUMPY_NULL
        self.assertAlmostEqual(np.abs(g_out-NUMPY_NULL).max(), 0.,
                msg="Wrong values in padding area")

    def test_numpy_data_2d(self):
        DIM=2
        testdata = np.arange(20*21).reshape(20,21)
        error = 1.*np.ones(testdata.shape)
        source = NumpyData(DataSource.GRAVITY, testdata, null_value=NUMPY_NULL)
        X0,NP,DX=source.getDataExtents()
        for i in range(DIM):
            self.assertAlmostEqual(X0[i], 0., msg="Data origin wrong")
            self.assertEqual(NP[i], testdata.shape[DIM-i-1], msg="Wrong number of data points")
            self.assertAlmostEqual(DX[i], 1000./testdata.shape[DIM-i-1], msg="Wrong cell size")

        domainbuilder=DomainBuilder(dim=3)
        domainbuilder.addSource(source)
        domainbuilder.setVerticalExtents(depth=-VMIN, air_layer=VMAX, num_cells=NE_V)
        domainbuilder.setElementPadding(PAD_X, PAD_Y)
        dom=domainbuilder.getDomain()
        g,s=domainbuilder.getGravitySurveys()[0]

        outfn=os.path.join(WORKDIR, '_npdata2d.csv')
        saveDataCSV(outfn, g=g, s=s)

        DV=(VMAX-VMIN)/NE_V

        # check data
        nx=NP[0]+2*PAD_X
        ny=NP[1]+2*PAD_Y
        nz=NE_V
        z_data=int(np.round((ALT-VMIN)/DV)-1)

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=float)
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)/(nx*ny)
        g_out=out[:,0].reshape(nz,ny,nx)
        s_out=out[:,1].reshape(nz,ny,nx)
        self.assertAlmostEqual(np.abs(
            g_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]-testdata).max(),
            0., msg="Difference in gravity data area")

        self.assertAlmostEqual(np.abs(
            s_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]-error).max(),
            0., msg="Difference in error data area")

        # overwrite data -> should only be padding value left
        g_out[z_data, PAD_Y:PAD_Y+NP[1], PAD_X:PAD_X+NP[0]]=NUMPY_NULL
        self.assertAlmostEqual(np.abs(g_out-NUMPY_NULL).max(), 0.,
                msg="Wrong values in padding area")


@unittest.skipIf(not haveProj, 'pyproj not available')
class TestErMapperData(unittest.TestCase):
    def test_ers_with_padding(self):
        source = ErMapperData(DataSource.GRAVITY, headerfile=ERS_DATA, 
                              altitude=ALT)
        domainbuilder=DomainBuilder()
        domainbuilder.addSource(source)
        domainbuilder.setVerticalExtents(depth=-VMIN, air_layer=VMAX, num_cells=NE_V)
        domainbuilder.setElementPadding(PAD_X,PAD_Y)
        dom=domainbuilder.getDomain()
        g,s=domainbuilder.getGravitySurveys()[0]

        outfn=os.path.join(WORKDIR, '_ersdata.csv')
        saveDataCSV(outfn, g=g, s=s)

        X0,NP,DX=source.getDataExtents()
        DV=(VMAX-VMIN)/NE_V

        # check metadata
        self.assertEqual(NP, ERS_SIZE, msg="Wrong number of data points")
        # this test only works if gdal is available
        try:
            import osgeo.osr
            for i in range(len(ERS_ORIGIN)):
                self.assertAlmostEqual(X0[i], ERS_ORIGIN[i], msg="Data origin wrong")
        except ImportError:
            print("Skipping test of data origin since gdal is not installed.")

        # check data
        nx=NP[0]+2*PAD_X
        ny=NP[1]+2*PAD_Y
        nz=NE_V
        z_data=int(np.round((ALT-VMIN)/DV)-1)

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

@unittest.skipIf('NetCdfData' not in dir(), 'netCDF not available')
class TestNetCdfData(unittest.TestCase):
    def test_cdf_with_padding(self):
        source = NetCdfData(DataSource.GRAVITY, NC_DATA, ALT)
        domainbuilder=DomainBuilder()
        domainbuilder.addSource(source)
        domainbuilder.setVerticalExtents(depth=-VMIN, air_layer=VMAX, num_cells=NE_V)
        domainbuilder.setElementPadding(PAD_X,PAD_Y)
        dom=domainbuilder.getDomain()
        g,s=domainbuilder.getGravitySurveys()[0]

        outfn=os.path.join(WORKDIR, '_ncdata.csv')
        saveDataCSV(outfn, g=g, s=s)

        X0,NP,DX=source.getDataExtents()
        DV=(VMAX-VMIN)/NE_V

        # check metadata
        self.assertEqual(NP, NC_SIZE, msg="Wrong number of data points")
        # this only works if gdal is available
        try:
            import osgeo.osr
            for i in range(len(NC_ORIGIN)):
                self.assertAlmostEqual(X0[i], NC_ORIGIN[i], msg="Data origin wrong")
        except ImportError:
            print("Skipping test of data origin since gdal is not installed.")

        # check data
        nx=NP[0]+2*PAD_X
        ny=NP[1]+2*PAD_Y
        nz=NE_V
        z_data=int(np.round((ALT-VMIN)/DV)-1)

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
    suite.addTest(unittest.makeSuite(TestNumpyData))
    suite.addTest(unittest.makeSuite(TestErMapperData))
    suite.addTest(unittest.makeSuite(TestNetCdfData))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

