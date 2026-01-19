
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


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os
import sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle, Brick, oxleycpp
# from test_objects import Test_Dump, Test_SetDataPointValue, Test_saveCSV, Test_TableInterpolation
from test_objects import Test_Dump, Test_SetDataPointValue, Test_saveCSV
from test_objects import Test_Domain, Test_Lazy

from test_shared import Test_Shared

from run_escriptOnOxley import Test_SharedOnOxley, Test_DomainOnOxley, \
                        Test_TableInterpolationOnOxley, Test_DataOpsOnOxley, \
                        Test_CSVOnOxley


def test_Rectangle(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    return m

def test_Brick(**kwargs):
    m = Brick(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    return m

try:
     OXLEY_WORKDIR=os.environ['OXLEY_WORKDIR']
except KeyError:
     OXLEY_WORKDIR='.'

NE=4 # number elements, must be even
mpiSize=getMPISizeWorld()
for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
    NX=x
    NY=mpiSize//x
    if NX*NY == mpiSize:
        break

for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
    NXb=x[0]
    NYb=x[1]
    NZb=mpiSize//(x[0]*x[1])
    if NXb*NYb*NZb == mpiSize:
        break

class Test_SharedOnMultiOxley(Test_SharedOnOxley):
    def setUp(self):
        self.domain=test_Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.tol=0.001

    def tearDown(self):
        del self.domain
        del self.tol

class Test_DomainOnMultiOxley(Test_DomainOnOxley):
    def setUp(self):
        self.boundary_tag_list = [1, 2, 10, 20]
        self.domain=test_Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.rdomain=test_Rectangle(n0=(NE+6)*NX-1, n1=(NE+6)*NY-1, l0=1., l1=1., d0=NX, d1=NY)

    def tearDown(self):
        del self.domain
        del self.boundary_tag_list

class Test_DataOpsOnMultiOxley(Test_DataOpsOnOxley):
    def setUp(self):
        self.domain=test_Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.domain_with_different_number_of_samples=test_Rectangle(n0=7*NE*NX-1, n1=3*NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.domain_with_different_number_of_data_points_per_sample=test_Rectangle(n0=7*NE*NX-1, n1=3*NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.domain_with_different_sample_ordering=test_Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.filename_base=OXLEY_WORKDIR
        self.mainfs=Function(self.domain)
        self.otherfs=Solution(self.domain)

    def tearDown(self):
        del self.domain
        del self.domain_with_different_number_of_samples
        del self.domain_with_different_number_of_data_points_per_sample
        del self.domain_with_different_sample_ordering
        del self.mainfs
        del self.otherfs

@unittest.skipIf(mpiSize > 1, "Multiresolution domains require single process")
@unittest.skip("Oxley multi-resolution Brick causes heap corruption - see issue #118")
class Test_TableInterpolationOnMultiOxley(Test_TableInterpolationOnOxley):
    def setUp(self):
        self.domain = test_Brick(n0=NE*NXb-1, n1=NE*NYb-1, n2=NE*NZb-1, l0=1., l1=1., l2=1., d0=NXb, d1=NYb, d2=NZb)
        self.functionspaces=[ContinuousFunction(self.domain), Function(self.domain), ReducedFunction(self.domain),
            FunctionOnBoundary(self.domain), ReducedFunctionOnBoundary(self.domain)]
        #We aren't testing DiracDeltaFunctions
        self.xn=5 # number of grids on x axis
        self.yn=5 # number of grids on y axis
        self.zn=5

    def tearDown(self):
        del self.domain
        del self.functionspaces

class Test_CSVOnMultiOxley(Test_CSVOnOxley):

    # Skip CSV tests - the expected line counts are calculated for regular meshes,
    # but multi-resolution meshes have different data point counts
    @unittest.skip("Test expectations not configured for multi-resolution mesh structure - see issue #118")
    def test_saveCSV_functionspaces(self):
        super().test_saveCSV_functionspaces()

    def setUp(self):
        self.workdir=OXLEY_WORKDIR
        self.domain=test_Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.functionspaces=[ContinuousFunction, Function, ReducedFunction,
                             FunctionOnBoundary, ReducedFunctionOnBoundary]

        NE0 = (NE*NX-1)*2
        NE1 = (NE*NY-1)*2

        # number of total data points for each function space
        self.linecounts=[ (NE0+1)*(NE1+1)+1, 4*NE0*NE1+1, NE0*NE1+1,
                4*NE0+4*NE1+1, 2*NE0+2*NE1+1 ]
        # number of masked points, i.e. where X[0] is non-zero
        self.linecounts_masked=[ NE0*(NE1+1)+1, 4*NE0*NE1+1, NE0*NE1+1,
                4*NE0+2*NE1+1, 2*NE0+NE1+1 ]
        # expected values in first line of masked data = [ X[:], X[0] ]
        self.firstline=[ [1./NE0, 0., 1./NE0],
                         [None, None, None],
                         [None, None, None],
                         [None, None, None],
                         [None, None, None] ]

    def tearDown(self):
        del self.domain


class Test_randomOnMultiOxley(unittest.TestCase):
    @unittest.skip("Oxley multi-resolution domains don't support Gaussian filter options - see issue #118")
    def test_FillRectangle(self):
        fs=ContinuousFunction(test_Rectangle(n0=5*(int(sqrt(mpiSize)+1)),n1=5*(int(sqrt(mpiSize)+1))))
        RandomData((), fs, 2,("gaussian",1,0.5))
        RandomData((), fs, 0,("gaussian",2,0.76))
        self.assertRaises(NotImplementedError, RandomData, (2,2), fs, 0, ("gaussian",2,0.76)) #data not scalar
        self.assertRaises(ValueError, RandomData, (), fs, 0, ("gaussian",11,0.1)) #radius too large
        RandomData((2,3),fs)

    @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
    @unittest.skip("Oxley multi-resolution Brick causes heap corruption - see issue #118")
    def test_FillBrick(self):
        # If we are going to do really big tests of this, the size of this brick will need to be reduced
        fs=ContinuousFunction(test_Brick(n0=5*mpiSize, n1=5*mpiSize, n2=5*mpiSize))
        RandomData((), fs, 2,("gaussian",1,0.5))
        RandomData((), fs, 0,("gaussian",2,0.76))
        self.assertRaises(NotImplementedError, RandomData, (2,2), fs, 0, ("gaussian",2,0.76)) #data not scalar
        self.assertRaises(ValueError, RandomData, (), fs, 0, ("gaussian",11,0.1)) #radius too large
        RandomData((2,3),fs)

# TODO
# class Test_multiResolution(unittest.TestCase):
#     # TODO
#     def test_MultiRectangle_constructors(self):
#         # with self.assertRaises(OverflowError): #negative is bad
#         #     MultiRectangle(n0=2*mpiSize-1, n1=5, d0=mpiSize, multiplier=-1)
#         # with self.assertRaises(RuntimeError): #zero is bad
#         #     MultiRectangle(n0=2*mpiSize-1, n1=5, d0=mpiSize, multiplier=0)
#         # with self.assertRaises(TypeError): #non-int is bad
#         #     MultiRectangle(n0=2*mpiSize-1, n1=5, d0=mpiSize, multiplier=.5)
#         # with self.assertRaises(RuntimeError): #non-power of two is bad
#         #     MultiRectangle(n0=2*mpiSize-1, n1=5, d0=mpiSize, multiplier=3)
#         # with self.assertRaises(Exception): #dimensions required
#         #     MultiRectangle(n1=5, d1=mpiSize, multiplier=3)
#         # MultiRectangle(n0=2*mpiSize-1, n1=5, d0=mpiSize, multiplier=1)
#         # MultiRectangle(n0=2*mpiSize-1, n1=5, d0=mpiSize, multiplier=2)
#         # MultiRectangle(n0=2*mpiSize-1, n1=5, d0=mpiSize, multiplier=8)
#         print("skipping...")

#     def test_RectangleInterpolation_NodesToNodesAndElements_CoarseToFine(self):
#         mrd = MultiResolutionDomain(2, n0=2, n1=2*mpiSize-1, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [i.getX() for i in domains]

#         for targetFS, name in [(Function, 'Function'),
#                        (ContinuousFunction, 'ContinuousFunction'),
#                        (ReducedContinuousFunction, 'ReducedContinuousFunction')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(source_level + 1, len(domains)):
#                     val = Lsup(interpolate(X[target_level], targetFS(domains[target_level])) \
#                             - interpolate(X[source_level], targetFS(domains[target_level])))
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'ContinuousFunction', source_level, name, target_level, val))

#     def test_RectangleInterpolation_NodesToElements_FineToCoarse(self):
#         mrd = MultiResolutionDomain(2, n0=2, n1=2*mpiSize-1, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [i.getX() for i in domains]

#         for targetFS, name in [(Function, 'Function')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(0, source_level):
#                     val = Lsup(interpolate(X[target_level], targetFS(domains[target_level])) \
#                             - interpolate(X[source_level], targetFS(domains[target_level])))
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'ContinuousFunction', source_level, name, target_level, val))

#     def test_RectangleInterpolation_ReducedToElements_CoarseToFine(self):
#         mrd = MultiResolutionDomain(2, n0=2, n1=2*mpiSize-1, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [interpolate(i.getX(), ReducedFunction(i)) for i in domains]

#         for targetFS, name in [(Function, 'Function')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(source_level + 1, len(domains)):
#                     to = targetFS(domains[target_level])
#                     desired = interpolate(X[source_level], Function(domains[source_level]))
#                     val = Lsup(interpolate(X[source_level], to) - desired)
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'ReducedFunction', source_level, name, target_level, val))        

#     def test_RectangleInterpolation_ElementsToElements_CoarseToFine(self):
#         mrd = MultiResolutionDomain(2, n0=2, n1=2*mpiSize-1, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [interpolate(i.getX(), Function(i)) for i in domains]

#         for targetFS, name in [(Function, 'Function')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(source_level + 1, len(domains)):
#                     val = Lsup(interpolate(X[source_level], targetFS(domains[target_level])) \
#                             - interpolate(X[target_level], targetFS(domains[target_level])))
#                     if val > 1e-12:
#                         print("Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'Function', source_level, name, target_level, val))
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'Function', source_level, name, target_level, val))

#     def test_RectangleInterpolation_ElementsToElements_FineToCoarse(self):
#         mrd = MultiResolutionDomain(2, n0=2, n1=2*mpiSize-1, d1=mpiSize, l0=2)
#         d0 = mrd.getLevel(0)
#         d1 = mrd.getLevel(1)
#         d2 = mrd.getLevel(2)
#         x0 = interpolate(d0.getX(), Function(d0))
#         x1 = interpolate(d1.getX(), Function(d1))
#         x2 = interpolate(d2.getX(), Function(d2))

#         val = Lsup(x0 - interpolate(x1, Function(d0)))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure from level 1 to level 0: %g !< 1e-12"%val)

#         val = Lsup(x1 - interpolate(x2, Function(d1)))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure from level 2 to level 1: %g !< 1e-12"%val)

#         val = Lsup(x0 - interpolate(x2, Function(d0)))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure from level 2 to level 0: %g !< 1e-12"%val)
        
#         val = Lsup(integrate(interpolate(sin(x2[0]), Function(d0))*x0) - integrate(sin(x2[0])*x2))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure: %g !< 1e-12"%val)

#         val = integrate(interpolate(sin(x2[0]), Function(d0))*x0[0]*x0[1]) - integrate(sin(x2[0])*x2[0]*x2[1])
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure: %g !< 1e-12"%val)
        
#         val = integrate(interpolate(sin(x2[0]), Function(d0))) - integrate(sin(x2[0]))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure: %g !< 1e-12"%val)



#     @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
#     def test_MultiBrick_constructors(self):
#         with self.assertRaises(OverflowError): #negative is bad
#             MultiBrick(n0=2*mpiSize-1, n1=5, n2=3, d1=mpiSize, multiplier=-1)
#         with self.assertRaises(RuntimeError): #zero is bad
#             MultiBrick(n0=2*mpiSize-1, n1=5, n2=3, d1=mpiSize, multiplier=0)
#         with self.assertRaises(TypeError): #non-int is bad
#             MultiBrick(n0=2*mpiSize-1, n1=5, n2=3, d1=mpiSize, multiplier=.5)
#         with self.assertRaises(RuntimeError): #non-power of two is bad
#             MultiBrick(n0=2*mpiSize-1, n1=5, n2=3, d1=mpiSize, multiplier=3)
#         with self.assertRaises(Exception): #dimensions required
#             MultiBrick(n1=5, n2=3, d1=mpiSize, multiplier=3)
#         MultiBrick(n0=2*mpiSize-1, n1=5, n2=3, d1=mpiSize, multiplier=1)
#         MultiBrick(n0=2*mpiSize-1, n1=5, n2=3, d1=mpiSize, multiplier=2)
#         MultiBrick(n0=2*mpiSize-1, n1=5, n2=3, d1=mpiSize, multiplier=8)

#     @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
#     def test_BrickInterpolation_NodesToNodesAndElements_CoarseToFine(self):
#         mrd = MultiResolutionDomain(3, n0=2, n1=2*mpiSize, n2=3, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [i.getX() for i in domains]

#         for targetFS, name in [(Function, 'Function'),
#                        (ContinuousFunction, 'ContinuousFunction'),
#                        (ReducedContinuousFunction, 'ReducedContinuousFunction')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(source_level + 1, len(domains)):
#                     val = Lsup(interpolate(X[target_level], targetFS(domains[target_level])) \
#                             - interpolate(X[source_level], targetFS(domains[target_level])))
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'ContinuousFunction', source_level, name, target_level, val))
#     @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
#     def test_BrickInterpolation_NodesToElements_FineToCoarse(self):
#         mrd = MultiResolutionDomain(3, n0=2, n1=2*mpiSize, n2=3, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [i.getX() for i in domains]

#         for targetFS, name in [(Function, 'Function')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(0, source_level):
#                     val = Lsup(interpolate(X[target_level], targetFS(domains[target_level])) \
#                             - interpolate(X[source_level], targetFS(domains[target_level])))
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'ContinuousFunction', source_level, name, target_level, val))
#     @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
#     def test_BrickInterpolation_ReducedToElements_CoarseToFine(self):
#         mrd = MultiResolutionDomain(3, n0=2, n1=2*mpiSize, n2=3, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [interpolate(i.getX(), ReducedFunction(i)) for i in domains]

#         for targetFS, name in [(Function, 'Function')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(source_level + 1, len(domains)):
#                     to = targetFS(domains[target_level])
#                     desired = interpolate(X[source_level], Function(domains[source_level]))
#                     val = Lsup(interpolate(X[source_level], to) \
#                             - interpolate(desired, to))
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'ReducedFunction', source_level, name, target_level, val))        

#     @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
#     def test_BrickInterpolation_ElementsToElements_CoarseToFine(self):
#         mrd = MultiResolutionDomain(3, n0=2, n1=2*mpiSize, n2=2, d1=mpiSize, l0=2)
#         domains = [mrd.getLevel(i) for i in range(3)]
#         X = [interpolate(i.getX(), Function(i)) for i in domains]

#         for targetFS, name in [(Function, 'Function')]:
#             for source_level in range(len(domains)):
#                 for target_level in range(source_level + 1, len(domains)):
#                     val = Lsup(interpolate(X[target_level], targetFS(domains[target_level])) \
#                             - interpolate(X[source_level], targetFS(domains[target_level])))
#                     self.assertLess(val, 1e-12,
#                             "Interpolation failure from %s level %d to %s level %d: %g !< 1e-12"%(\
#                             'Function', source_level, name, target_level, val))

#     @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
#     def test_BrickInterpolation_ElementsToElements_FineToCoarse(self):
#         mrd = MultiResolutionDomain(3, n0=2, n1=2*mpiSize, n2=3, d1=mpiSize, l0=2)
#         d0 = mrd.getLevel(0)
#         d1 = mrd.getLevel(1)
#         d2 = mrd.getLevel(2)
#         x0 = interpolate(d0.getX(), Function(d0))
#         x1 = interpolate(d1.getX(), Function(d1))
#         x2 = interpolate(d2.getX(), Function(d2))

#         val = Lsup(x0 - interpolate(x1, Function(d0)))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure from level 1 to level 0: %g !< 1e-12"%val)

#         val = Lsup(x1 - interpolate(x2, Function(d1)))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure from level 2 to level 1: %g !< 1e-12"%val)

#         val = Lsup(x0 - interpolate(x2, Function(d0)))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure from level 2 to level 0: %g !< 1e-12"%val)
        
#         val = Lsup(integrate(interpolate(sin(x2[0]), Function(d0))*x0) - integrate(sin(x2[0])*x2))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure: %g !< 1e-12"%val)

#         val = integrate(interpolate(sin(x2[0]), Function(d0))*x0[0]*x0[1]*x0[2]) - integrate(sin(x2[0])*x2[0]*x2[1]*x2[2])
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure: %g !< 1e-12"%val)
        
#         val = integrate(interpolate(sin(x2[0]), Function(d0))) - integrate(sin(x2[0]))
#         self.assertLess(val, 1e-12,
#                 "Interpolation failure: %g !< 1e-12"%val)



if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

