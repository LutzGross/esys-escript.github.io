
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

import os
import numpy as np
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.ripley import Rectangle, Brick, ripleycpp

try:
     RIPLEY_WORKDIR=os.environ['RIPLEY_WORKDIR']
except KeyError:
     RIPLEY_WORKDIR='.'

#NE=4 # number elements, must be even
#for x in [int(sqrt(ranks)),2,3,5,7,1]:
#    NX=x
#    NY=ranks//x
#    if NX*NY == ranks:
#        break
#
#for x in [(int(ranks**(1/3.)),int(ranks**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
#    NXb=x[0]
#    NYb=x[1]
#    NZb=ranks//(x[0]*x[1])
#    if NXb*NYb*NZb == ranks:
#        break

ranks = getMPISizeWorld()

def adjust(NE, ftype):
    if ftype == ContinuousFunction:
        return [i+1 for i in NE]
    return NE


class WriteBinaryGridTestBase(unittest.TestCase):
    NX = 10*ranks-1
    NZ = 10

    def write(self, domain, data, filename):
        domain.writeBinaryGrid(data, filename, self.byteorder, self.datatype)

    def generateUniqueData(self, ftype):
        dim = self.domain.getDim()
        FSx=ftype(self.domain).getX()
        NE = adjust(self.NE, ftype)
        # normalise and scale range of values
        x = [FSx[i]-inf(FSx[i]) for i in range(dim)]
        x = [(NE[i]-1)*(x[i]/sup(x[i])) for i in range(dim)]
        xMax = [int(sup(x[i]))+1 for i in range(dim)]
        nvals=NE[0]*NE[1]
        data = x[0] + xMax[0]*x[1]
        if self.datatype == ripleycpp.DATATYPE_INT32:
            data += 0.05
        if dim > 2:
            data = data + xMax[0]*xMax[1]*x[2]
            nvals*=NE[2]

        grid = np.array(range(nvals), dtype=self.dtype).reshape(tuple(reversed(NE)))
        return data, grid

    def test_writeGrid2D(self):
        self.NE = [self.NX, self.NZ]
        self.domain = Rectangle(self.NE[0], self.NE[1], d1=0)
        for ftype in [ReducedFunction, ContinuousFunction]:
            data, original_grid = self.generateUniqueData(ftype)
            filename = os.path.join(RIPLEY_WORKDIR, "tempfile")
            self.write(self.domain, data, filename)
            MPIBarrierWorld()
            result = np.fromfile(filename, dtype=self.dtype).reshape(
                    tuple(reversed(adjust(self.NE,ftype))))
            self.assertEquals(Lsup(original_grid-result), 0, "Data doesn't match for "+str(ftype(self.domain)))

    def test_writeGrid3D(self):
        self.NE = [self.NX, self.NX, self.NZ]
        self.domain = Brick(self.NE[0], self.NE[1], self.NE[2], d2=0)
        for ftype in [ReducedFunction, ContinuousFunction]:
            data, original_grid = self.generateUniqueData(ftype)
            filename = os.path.join(RIPLEY_WORKDIR, "tempfile")
            self.write(self.domain, data, filename)
            MPIBarrierWorld()
            result = np.fromfile(filename, dtype=self.dtype).reshape(
                    tuple(reversed(adjust(self.NE,ftype))))
            self.assertEquals(Lsup(original_grid-result), 0, "Data doesn't match for "+str(ftype(self.domain)))

class Test_writeBinaryGridRipley_LITTLE_FLOAT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = ripleycpp.BYTEORDER_LITTLE_ENDIAN
        self.datatype = ripleycpp.DATATYPE_FLOAT32
        self.dtype = "<f4"

class Test_writeBinaryGridRipley_LITTLE_FLOAT64(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = ripleycpp.BYTEORDER_LITTLE_ENDIAN
        self.datatype = ripleycpp.DATATYPE_FLOAT64
        self.dtype = "<f8"

class Test_writeBinaryGridRipley_LITTLE_INT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = ripleycpp.BYTEORDER_LITTLE_ENDIAN
        self.datatype = ripleycpp.DATATYPE_INT32
        self.dtype = "<i4"

class Test_writeBinaryGridRipley_BIG_FLOAT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = ripleycpp.BYTEORDER_BIG_ENDIAN
        self.datatype = ripleycpp.DATATYPE_FLOAT32
        self.dtype = ">f4"

class Test_writeBinaryGridRipley_BIG_FLOAT64(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = ripleycpp.BYTEORDER_BIG_ENDIAN
        self.datatype = ripleycpp.DATATYPE_FLOAT64
        self.dtype = ">f8"

class Test_writeBinaryGridRipley_BIG_INT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = ripleycpp.BYTEORDER_BIG_ENDIAN
        self.datatype = ripleycpp.DATATYPE_INT32
        self.dtype = ">i4"


@unittest.skipIf(getMPISizeWorld() > 1,
    "Skipping compressed binary grid tests due to element stretching")
class Test_readBinaryGridZippedRipley(unittest.TestCase):
    # constants
    byteorder = ripleycpp.BYTEORDER_NATIVE
    datatype = ripleycpp.DATATYPE_FLOAT64

    def read(self, filename, FS, expected, zipped = False):
        first = [0 for i in expected]
        reverse = [0 for i in expected]
        scale = [1 for i in expected]

        if not zipped:
            return ripleycpp._readBinaryGrid(filename, FS, (), 50000,
                self.byteorder, self.datatype, first, expected, scale, reverse)

        if not hasattr(ripleycpp, "_readBinaryGridFromZipped"):
            raise unittest.SkipTest("unzip library not available (boost_iostreams)")
        return ripleycpp._readBinaryGridFromZipped(filename, FS, (), 50000,
                self.byteorder, self.datatype, first, expected, scale, reverse)

    def test_readCompressed2D(self):
        NE = [9, 10]
        domain = Rectangle(NE[0], NE[1], d1=0)
        for filename, ftype in [("RectRedF%s.grid.gz", ReducedFunction),
                ("RectConF%s.grid.gz", ContinuousFunction)]:
            FS = ftype(domain)
            filename = os.path.join("ref_data", filename%ranks)
            print(filename)
            unzipped = self.read(filename[:-3], FS, adjust(NE, ftype))
            zipped = self.read(filename, FS, adjust(NE, ftype), True)
            self.assertEqual(Lsup(zipped - unzipped), 0, "Data objects don't match for "+str(FS))

    def test_readCompressed3D(self):
        NE = [9, 9, 10]
        domain = Brick(NE[0], NE[1], NE[2], d1=0, d2=0)
        for filename, ftype in [("BrickRedF%s.grid.gz", ReducedFunction),
                ("BrickConF%s.grid.gz", ContinuousFunction)]:
            FS = ftype(domain)
            filename = os.path.join("ref_data", filename%ranks)
            unzipped = self.read(filename[:-3], FS, adjust(NE, ftype))
            zipped = self.read(filename, FS, adjust(NE, ftype), True)
            self.assertEqual(Lsup(zipped - unzipped), 0, "Data objects don't match for "+str(FS))


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

