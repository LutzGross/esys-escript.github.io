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
import numpy as np
import os
import sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import inf, sup, saveDataCSV, getMPISizeWorld, hasFeature
from esys.downunder.datasources import *
from esys.downunder.domainbuilder import DomainBuilder
from esys.downunder.coordinates import WGS84ReferenceSystem

HAVE_RIPLEY = True
try:
    from esys.ripley import Rectangle
except ImportError as e:
    HAVE_RIPLEY = False

if hasFeature('netcdf'):
    HAVE_NETCDF = True
else:
    HAVE_NETCDF = False

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

ERS32_DATA = os.path.join(TEST_DATA_ROOT, 'ermapper32_test.ers')
ERS64_DATA = os.path.join(TEST_DATA_ROOT, 'ermapper64_test.ers')
ERS_REF = os.path.join(TEST_DATA_ROOT, 'ermapper_test.csv')
ERS_NULL = -99999 * 1e-6
ERS_SIZE = [20,15]
ERS_ORIGIN = [309241.0, 6318655.0]
NC_DATA = os.path.join(TEST_DATA_ROOT, 'netcdf_test.nc')
NC_REF = os.path.join(TEST_DATA_ROOT, 'netcdf_test.csv')
NC_NULL = 0.
NC_SIZE = [20,15]
NC_ORIGIN = [403320.91466610413, 6414860.942530109]
NC_ORIGIN_WGS84 = [115.97200299999747, -32.399100000001042]
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
    @unittest.skipIf(not HAVE_RIPLEY, "ripley module not available")
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

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=np.float64)
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)//nx
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

    @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
    @unittest.skipIf(not HAVE_RIPLEY, "ripley module not available")
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

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=np.float64)
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)//(nx*ny)
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
@unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
class TestErMapperData(unittest.TestCase):
    @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
    def test_ers32_with_padding(self):
        self._ers_tester(ERS32_DATA)

    @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
    def test_ers64_with_padding(self):
        self._ers_tester(ERS64_DATA)

    def _ers_tester(self, filename):
        source = ErMapperData(DataSource.GRAVITY, headerfile=filename, 
                              altitude=ALT, scale_factor=1e-6)
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

        ref=np.genfromtxt(ERS_REF, delimiter=',', dtype=np.float64)
        g_ref=ref[:,0].reshape((NP[1],NP[0]))
        s_ref=ref[:,1].reshape((NP[1],NP[0]))

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=np.float64)
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)//(nx*ny)
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

@unittest.skipIf(HAVE_NETCDF is False, "netCDF not available")
@unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
@unittest.skipIf('NetCdfData' not in dir(), 'netCDF not available')
class TestNetCdfData(unittest.TestCase):
    @unittest.skipIf(not haveProj, 'pyproj not available')
    @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
    def test_cdf_with_padding(self):
        source = NetCdfData(DataSource.GRAVITY, NC_DATA, ALT, scale_factor=1e-6)
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
                self.assertAlmostEqual(X0[i], NC_ORIGIN[i], places=3, msg="Data origin wrong")
        except ImportError:
            print("Skipping test of data origin since gdal is not installed.")

        # check data
        nx=NP[0]+2*PAD_X
        ny=NP[1]+2*PAD_Y
        nz=NE_V
        z_data=int(np.round((ALT-VMIN)/DV)-1)
    
        ref=np.genfromtxt(NC_REF, delimiter=',', dtype=np.float64)
        g_ref=ref[:,0].reshape((NP[1],NP[0]))
        s_ref=ref[:,1].reshape((NP[1],NP[0]))

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=np.float64)
        # recompute nz since ripley might have adjusted number of elements

        nz=len(out)//(nx*ny)
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

    @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
    def test_cdf_with_padding_ellipsoid(self):
        ref=WGS84ReferenceSystem()

        source = NetCdfData(DataSource.GRAVITY, NC_DATA, ALT,
                            reference_system=ref, scale_factor=1e-6)
        domainbuilder=DomainBuilder(reference_system=ref)
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

        for i in range(len(NC_ORIGIN)):
            self.assertAlmostEqual(X0[i], NC_ORIGIN_WGS84[i], msg="Data origin wrong")

        # check data
        nx=NP[0]+2*PAD_X
        ny=NP[1]+2*PAD_Y
        nz=NE_V
        z_data=int(np.round((ALT-VMIN)/DV)-1)

        ref=np.genfromtxt(NC_REF, delimiter=',', dtype=np.float64)
        g_ref=ref[:,0].reshape((NP[1],NP[0]))
        s_ref=ref[:,1].reshape((NP[1],NP[0]))

        out=np.genfromtxt(outfn, delimiter=',', skip_header=1, dtype=np.float64)
        
        # recompute nz since ripley might have adjusted number of elements
        nz=len(out)//(nx*ny)
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

class TestSeismicSource(unittest.TestCase):
    def test_seismic_source(self):
        ss= SeismicSource(x=1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,6.], elevation=7.)
        self.assertEqual(ss.getLocation(), (1.,2) )
        self.assertEqual(ss.getFrequency(), 3.)
        self.assertEqual(ss.getPower(), complex(4.,-4.) ) 
        self.assertEqual(ss.getOrientation(), [5.,6.])
        self.assertEqual(ss.getElevation(), 7.)
        
        self.assertTrue( ss == SeismicSource(x=1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertFalse( ss == SeismicSource(x=-1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertFalse( ss == SeismicSource(x=1., y=-2., omega=3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertFalse( ss == SeismicSource(x=1., y=2., omega=-3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertFalse( ss == SeismicSource(x=1., y=2., omega=3., power = complex(-4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertFalse( ss == SeismicSource(x=1., y=2., omega=3., power = complex(4.,4.), orientation=[5.,6.],elevation=7.) )
        self.assertFalse( ss == SeismicSource(x=1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,-6.],elevation=7.) )
        self.assertFalse( ss == SeismicSource(x=1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,-6.],elevation=-7.) )

        self.assertFalse( ss != SeismicSource(x=1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertTrue( ss != SeismicSource(x=-1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertTrue( ss != SeismicSource(x=1., y=-2., omega=3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertTrue( ss != SeismicSource(x=1., y=2., omega=-3., power = complex(4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertTrue( ss != SeismicSource(x=1., y=2., omega=3., power = complex(-4.,-4.), orientation=[5.,6.],elevation=7.) )
        self.assertTrue( ss != SeismicSource(x=1., y=2., omega=3., power = complex(4.,4.), orientation=[5.,6.],elevation=7.) )
        self.assertTrue( ss != SeismicSource(x=1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,-6.],elevation=7.) )
        self.assertTrue( ss != SeismicSource(x=1., y=2., omega=3., power = complex(4.,-4.), orientation=[5.,-6.],elevation=-7.) )
        
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

