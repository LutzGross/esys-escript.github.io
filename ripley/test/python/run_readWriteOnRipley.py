
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import os
import numpy as np
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.ripley import *

try:
     RIPLEY_WORKDIR=os.environ['RIPLEY_WORKDIR']
except KeyError:
     RIPLEY_WORKDIR='/tmp'

HAVE_UNZIP = hasFeature('unzip')

#NE=4 # number elements, must be even
#for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
#    NX=x
#    NY=mpiSize//x
#    if NX*NY == mpiSize:
#        break
#
#for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
#    NXb=x[0]
#    NYb=x[1]
#    NZb=mpiSize//(x[0]*x[1])
#    if NXb*NYb*NZb == mpiSize:
#        break

mpiSize = getMPISizeWorld()
mpiRank = getMPIRankWorld()

def adjust(NE, ftype):
    if ftype in (ContinuousFunction, Solution):
        return [i+1 for i in NE]
    return NE


class WriteBinaryGridTestBase(unittest.TestCase): #subclassing required
    NX = 10*mpiSize-1
    NZ = 10

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
        if self.datatype == DATATYPE_INT32:
            data += 0.05
        if dim > 2:
            data = data + xMax[0]*xMax[1]*x[2]
            nvals*=NE[2]

        grid = np.array(range(nvals), dtype=self.dtype).reshape(tuple(reversed(NE)))
        return data, grid

    def writeThenRead(self, data, ftype, fcode):
        filename = os.path.join(RIPLEY_WORKDIR, "_wgrid%dd%s"%(self.domain.getDim(),fcode))
        filename = filename + self.dtype.replace('<','L').replace('>','B')
        self.domain.writeBinaryGrid(data, filename, self.byteorder, self.datatype)
        MPIBarrierWorld()
        result = np.fromfile(filename, dtype=self.dtype).reshape(
                tuple(reversed(adjust(self.NE,ftype))))
        return result

    def test_writeGrid2D(self):
        self.NE = [self.NX, self.NZ]
        self.domain = Rectangle(self.NE[0], self.NE[1], d1=0)
        for ftype,fcode in [(ReducedFunction,'RF'), (ContinuousFunction,'CF'), (Solution, 'Sol')]:
            data, ref = self.generateUniqueData(ftype)
            result = self.writeThenRead(data, ftype, fcode)
            self.assertAlmostEqual(Lsup(ref-result), 0, delta=1e-9,
                    msg="Data doesn't match for "+str(ftype(self.domain)))

    def test_writeGrid3D(self):
        self.NE = [self.NX, self.NX, self.NZ]
        self.domain = Brick(self.NE[0], self.NE[1], self.NE[2], d2=0)
        for ftype,fcode in [(ReducedFunction,'RF'), (ContinuousFunction,'CF'), (Solution, 'Sol')]:
            data, ref = self.generateUniqueData(ftype)
            result = self.writeThenRead(data, ftype, fcode)
            self.assertAlmostEqual(Lsup(ref-result), 0, delta=1e-9,
                    msg="Data doesn't match for "+str(ftype(self.domain)))

class Test_writeBinaryGridRipley_LITTLE_FLOAT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_LITTLE_ENDIAN
        self.datatype = DATATYPE_FLOAT32
        self.dtype = "<f4"

class Test_writeBinaryGridRipley_LITTLE_FLOAT64(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_LITTLE_ENDIAN
        self.datatype = DATATYPE_FLOAT64
        self.dtype = "<f8"

class Test_writeBinaryGridRipley_LITTLE_INT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_LITTLE_ENDIAN
        self.datatype = DATATYPE_INT32
        self.dtype = "<i4"

class Test_writeBinaryGridRipley_BIG_FLOAT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_BIG_ENDIAN
        self.datatype = DATATYPE_FLOAT32
        self.dtype = ">f4"

class Test_writeBinaryGridRipley_BIG_FLOAT64(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_BIG_ENDIAN
        self.datatype = DATATYPE_FLOAT64
        self.dtype = ">f8"

class Test_writeBinaryGridRipley_BIG_INT32(WriteBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_BIG_ENDIAN
        self.datatype = DATATYPE_INT32
        self.dtype = ">i4"


class ReadBinaryGridTestBase(unittest.TestCase): #subclassing required
    """
    The reader tests work in several stages:
    1) create numpy array and write to temporary file (ref)
    2) call readBinaryGrid with that filename
    3) write the resulting Data object using writeBinaryGrid (test)
    4) read the result using numpy and compare (ref) and (test)
    As such, it is important to note that a working writeBinaryGrid() method
    is assumed!
    """
    # set defaults which may be overridden in subclasses
    NX = 10
    NZ = 8
    fspaces = [(ReducedFunction,'RF'), (ContinuousFunction,'CF')]
    byteorder = BYTEORDER_NATIVE
    datatype = DATATYPE_FLOAT64
    dtype = "f8"
    shape = ()
    fill = -42.57
    first = [0,0,0]
    multiplier = [1,1,1]
    reverse = [0,0,0]

    def generateUniqueData(self, ftype):
        dim = self.domain.getDim()
        NE = adjust(self.Ndata, ftype)
        nvals=NE[0]*NE[1]
        if dim > 2:
            nvals*=NE[2]
        grid = np.array(range(nvals), dtype=self.dtype).reshape(tuple(reversed(NE)))
        return grid

    def write(self, data, filename):
        self.domain.writeBinaryGrid(data, filename, self.byteorder, self.datatype)

    def read(self, filename, ftype):
        first = self.first[:self.domain.getDim()]
        multiplier = self.multiplier[:self.domain.getDim()]
        reverse = self.reverse[:self.domain.getDim()]
        numValues=adjust(self.Ndata, ftype)
        return readBinaryGrid(filename, ftype(self.domain),
                shape=self.shape, fill=self.fill, byteOrder=self.byteorder,
                dataType=self.datatype, first=first, numValues=numValues,
                multiplier=multiplier, reverse=reverse)

    def numpy2Data2Numpy(self, ref, ftype, fcode):
        filename = os.path.join(RIPLEY_WORKDIR, "_rgrid%dd%s"%(self.domain.getDim(),fcode))
        filename = filename + self.dtype.replace('<','L').replace('>','B')
        if mpiRank == 0:
            ref.tofile(filename)
        MPIBarrierWorld()
        # step 2 - read
        data = self.read(filename, ftype)
        MPIBarrierWorld()
        # step 3 - write
        self.write(data, filename) # overwrite is ok
        MPIBarrierWorld()
        result = np.fromfile(filename, dtype=self.dtype).reshape(
                tuple(reversed(adjust(self.NE,ftype))))
        return result

    def test_readGrid2D(self):
        if self.multiplier[0] == 1:
            self.NE = [self.NX*mpiSize-1, self.NZ*self.multiplier[1]]
        else:
            self.NE = [self.NX*mpiSize*self.multiplier[0]-1, self.NZ*self.multiplier[1]]
        self.domain = Rectangle(self.NE[0], self.NE[1], d0=mpiSize, d1=1)
        for ftype,fcode in self.fspaces:
            self.Ndata = [self.NX*mpiSize-1, self.NZ]
            if ftype==ContinuousFunction:
                self.Ndata[1] = self.NZ-1
            # step 1 - generate
            ref = self.generateUniqueData(ftype)
            # step 2 & 3
            result = self.numpy2Data2Numpy(ref, ftype, fcode)
            # apply transformations to be able to compare
            if self.reverse[0]:
                result = result[...,::-1]
            if self.reverse[1]:
                result = result[::-1,:]
            for i in range(2):
                ref = np.repeat(ref, self.multiplier[i], axis=1-i)

            # if domain larger than data: add column(s)/row(s) with fill value
            fill=np.array(self.fill, dtype=ref.dtype)
            realNE = adjust(self.NE,ftype)
            for d in range(2):
                excess = realNE[d]-ref.shape[1-d]
                if excess > 0:
                    shape = list(ref.shape)
                    shape[1-d] = excess
                    extra = fill * np.ones(shape)
                    if self.reverse[d]:
                        ref = np.append(extra, ref, axis=1-d)
                    else:
                        ref = np.append(ref, extra, axis=1-d)

            # step 4 - compare
            self.assertAlmostEqual(Lsup(ref-result), 0, delta=1e-9,
                    msg="Data doesn't match for "+str(ftype(self.domain)))

    def test_readGrid3D(self):
        if self.multiplier[0] == 1:
            self.NE = [self.NX*mpiSize-1, self.NX*self.multiplier[1], self.NZ*self.multiplier[2]]
        else:
            self.NE = [self.NX*mpiSize*self.multiplier[0]-1,
                       self.NX*self.multiplier[1], self.NZ*self.multiplier[2]]
        self.domain = Brick(self.NE[0], self.NE[1], self.NE[2], d0=mpiSize, d1=1, d2=1)
        for ftype,fcode in self.fspaces:
            self.Ndata = [self.NX*mpiSize-1, self.NX, self.NZ]
            if ftype==ContinuousFunction:
                self.Ndata[1] = self.NX-1
                self.Ndata[2] = self.NZ-1
            # step 1 - generate
            ref = self.generateUniqueData(ftype)
            # step 2 & 3
            result = self.numpy2Data2Numpy(ref, ftype, fcode)
            # apply transformations to be able to compare
            if self.reverse[0]:
                result = result[...,::-1]
            if self.reverse[1]:
                result = result[...,::-1,:]
            if self.reverse[2]:
                result = result[::-1,:,:]
            for i in range(3):
                ref = np.repeat(ref, self.multiplier[i], axis=2-i)

            # if domain larger than data: add column(s)/row(s) with fill value
            fill=np.array(self.fill, dtype=ref.dtype)
            realNE = adjust(self.NE,ftype)
            for d in range(3):
                excess = realNE[d]-ref.shape[2-d]
                if excess > 0:
                    shape = list(ref.shape)
                    shape[2-d] = excess
                    extra = fill * np.ones(shape)
                    if self.reverse[d]:
                        ref = np.append(extra, ref, axis=2-d)
                    else:
                        ref = np.append(ref, extra, axis=2-d)

            # step 4 - compare
            self.assertAlmostEqual(Lsup(ref-result), 0, delta=1e-9,
                    msg="Data doesn't match for "+str(ftype(self.domain)))


# The following block tests the reader for different byte orders and data
# types with domain-filling data (i.e. multiplier=1, reverse=0 and N=NE)

class Test_readBinaryGridRipley_LITTLE_FLOAT32(ReadBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_LITTLE_ENDIAN
        self.datatype = DATATYPE_FLOAT32
        self.dtype = "<f4"

class Test_readBinaryGridRipley_LITTLE_FLOAT64(ReadBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_LITTLE_ENDIAN
        self.datatype = DATATYPE_FLOAT64
        self.dtype = "<f8"

class Test_readBinaryGridRipley_LITTLE_INT32(ReadBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_LITTLE_ENDIAN
        self.datatype = DATATYPE_INT32
        self.dtype = "<i4"

class Test_readBinaryGridRipley_BIG_FLOAT32(ReadBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_BIG_ENDIAN
        self.datatype = DATATYPE_FLOAT32
        self.dtype = ">f4"

class Test_readBinaryGridRipley_BIG_FLOAT64(ReadBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_BIG_ENDIAN
        self.datatype = DATATYPE_FLOAT64
        self.dtype = ">f8"

class Test_readBinaryGridRipley_BIG_INT32(ReadBinaryGridTestBase):
    def setUp(self):
        self.byteorder = BYTEORDER_BIG_ENDIAN
        self.datatype = DATATYPE_INT32
        self.dtype = ">i4"

@unittest.skip("reverseX not supported yet")
class Test_readBinaryGridRipley_reverseX(ReadBinaryGridTestBase):
    def setUp(self):
        self.reverse = [1,0,0]

@unittest.skip("reverseY not supported yet")
class Test_readBinaryGridRipley_reverseY(ReadBinaryGridTestBase):
    def setUp(self):
        self.reverse = [0,1,0]

class Test_readBinaryGridRipley_reverseZ(ReadBinaryGridTestBase):
    def setUp(self):
        self.reverse = [0,0,1]

class Test_readBinaryGridRipley_multiplierX(ReadBinaryGridTestBase):
    def setUp(self):
        self.multiplier = [2,1,1]

class Test_readBinaryGridRipley_multiplierY(ReadBinaryGridTestBase):
    def setUp(self):
        self.multiplier = [1,2,1]

class Test_readBinaryGridRipley_multiplierZ(ReadBinaryGridTestBase):
    def setUp(self):
        self.multiplier = [1,1,2]

class Test_readBinaryGridRipley_multiplierXYZ(ReadBinaryGridTestBase):
    def setUp(self):
        self.multiplier = [2,3,4]


@unittest.skipIf(getMPISizeWorld() > 1,
    "Skipping compressed binary grid tests due to element stretching")
class Test_readBinaryGridZippedRipley(unittest.TestCase):
    # constants
    byteorder = BYTEORDER_NATIVE
    datatype = DATATYPE_FLOAT64

    def read(self, filename, FS, expected, zipped = False):
        first = [0 for i in expected]
        reverse = [0 for i in expected]
        scale = [1 for i in expected]

        if not zipped:
            return readBinaryGrid(filename, FS, (), 50000,
                self.byteorder, self.datatype, first, expected, scale, reverse)

        if not HAVE_UNZIP:
            raise unittest.SkipTest("unzip library not available (boost_iostreams)")
        return ripleycpp._readBinaryGridFromZipped(filename, FS, (), 50000,
                self.byteorder, self.datatype, first, expected, scale, reverse)

    def test_readCompressed2D(self):
        NE = [9, 10]
        domain = Rectangle(NE[0], NE[1], d1=0)
        for filename, ftype in [("RectRedF%s.grid.gz", ReducedFunction),
                ("RectConF%s.grid.gz", ContinuousFunction)]:
            FS = ftype(domain)
            filename = os.path.join("ref_data", filename%mpiSize)
            unzipped = self.read(filename[:-3], FS, adjust(NE, ftype))
            zipped = self.read(filename, FS, adjust(NE, ftype), True)
            self.assertEqual(Lsup(zipped - unzipped), 0, "Data objects don't match for "+str(FS))

    def test_readCompressed3D(self):
        NE = [9, 9, 10]
        domain = Brick(NE[0], NE[1], NE[2], d1=0, d2=0)
        for filename, ftype in [("BrickRedF%s.grid.gz", ReducedFunction),
                ("BrickConF%s.grid.gz", ContinuousFunction)]:
            FS = ftype(domain)
            filename = os.path.join("ref_data", filename%mpiSize)
            unzipped = self.read(filename[:-3], FS, adjust(NE, ftype))
            zipped = self.read(filename, FS, adjust(NE, ftype), True)
            self.assertEqual(Lsup(zipped - unzipped), 0, "Data objects don't match for "+str(FS))


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

