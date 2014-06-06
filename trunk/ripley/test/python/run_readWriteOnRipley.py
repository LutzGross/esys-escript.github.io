
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
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.ripley import Rectangle, Brick, ripleycpp

try:
     RIPLEY_WORKDIR=os.environ['RIPLEY_WORKDIR']
except KeyError:
     RIPLEY_WORKDIR='.'

#NE=4 # number elements, must be even
#mpiSize=getMPISizeWorld()
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

class Test_writeBinaryGridOnRipley(unittest.TestCase):
    def setUp(self):
        self.byteorder = ripleycpp.BYTEORDER_NATIVE
        self.datatype = ripleycpp.DATATYPE_FLOAT64
        self.ranks = getMPISizeWorld()

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

    def write(self, domain, data, filename):
        domain.writeBinaryGrid(data, filename, self.byteorder, self.datatype)

    def adjust(self, NE, ftype):
        if ftype == ContinuousFunction:
            return [i+1 for i in NE]
        return NE

    def generateUniqueData(self, FS, dim):
        x = FS.getX()
        if dim == 2:
            return x[0] * 10 * (10*self.ranks-1) + x[1]
        return x[0] * 100 * (10*self.ranks-1) + x[1] * 100 + x[2]

    def test_BrickWriteThenRead(self):
        NE = [10*self.ranks-1, 10*self.ranks-1, 10]
        domain = Brick(NE[0], NE[1], NE[2], d2=0)
        for ftype in [ReducedFunction, ContinuousFunction]:
            FS = ftype(domain)
            original = self.generateUniqueData(FS, 3)
            filename = os.path.join(RIPLEY_WORKDIR, "tempfile")
            self.write(domain, original, filename)
            result = self.read(filename, FS, self.adjust(NE, ftype))
            MPIBarrierWorld()
            if getMPIRankWorld() == 0:
                os.unlink(filename)
            self.assertEqual(Lsup(original - result), 0, "Data objects don't match for "+str(FS))

    def test_RectangleWriteThenRead(self):
        NE = [10*self.ranks-1, 10]
        domain = Rectangle(NE[0], NE[1], d1=0)
        for ftype in [ReducedFunction, ContinuousFunction]:
            FS = ftype(domain)
            original = self.generateUniqueData(FS, 2)
            filename = os.path.join(RIPLEY_WORKDIR, "tempfile")
            self.write(domain, original, filename)
            result = self.read(filename, FS, self.adjust(NE, ftype))
            MPIBarrierWorld()
            if getMPIRankWorld() == 0:
                os.unlink(filename)
            self.assertEqual(Lsup(original - result), 0, "Data objects don't match for "+str(FS))

    @unittest.skipIf(getMPISizeWorld() > 1,
        "Skipping compressed binary grid tests due to element stretching")
    def test_BrickCompressed(self):
        NE = [10*self.ranks-1, 10*self.ranks-1, 10]
        domain = Brick(NE[0], NE[1], NE[2], d1=0, d2=0)
        for filename, ftype in [("BrickRedF%s.grid.gz", ReducedFunction),
                ("BrickConF%s.grid.gz", ContinuousFunction)]:
            FS = ftype(domain)
            filename = os.path.join('ref_data', filename%self.ranks)
            unzipped = self.read(filename[:-3], FS, self.adjust(NE, ftype))
            zipped = self.read(filename, FS, self.adjust(NE, ftype), zipped=True)
            self.assertEqual(Lsup(zipped - unzipped), 0, "Data objects don't match for "+str(FS))

    @unittest.skipIf(getMPISizeWorld() > 1,
        "Skipping compressed binary grid tests due element stretching")
    def test_RectangleCompressed(self):
        NE = [10*self.ranks-1, 10]
        domain = Rectangle(NE[0], NE[1], d1=0)
        for filename, ftype in [("ref_data/RectRedF%s.grid.gz", ReducedFunction),
                ("ref_data/RectConF%s.grid.gz", ContinuousFunction)]:
            FS = ftype(domain)
            filename = filename%self.ranks
            unzipped = self.read(filename[:-3], FS, self.adjust(NE, ftype))
            zipped = self.read(filename, FS, self.adjust(NE, ftype), zipped=True)
            self.assertEqual(Lsup(zipped - unzipped), 0, "Data objects don't match for "+str(FS))

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

