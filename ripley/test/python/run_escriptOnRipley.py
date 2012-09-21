
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

import unittest
import tempfile

from esys.escript import *
from esys.ripley import Rectangle, Brick
import sys
import os
from test_objects import Test_Dump, Test_SetDataPointValue, Test_saveCSV, Test_TableInterpolation
from test_objects import Test_Domain, Test_GlobalMinMax, Test_Lazy

from test_shared import Test_Shared

try:
     RIPLEY_WORKDIR=os.environ['RIPLEY_WORKDIR']
except KeyError:
     RIPLEY_WORKDIR='.'

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

class Test_SharedOnRipley(Test_Shared):
    def setUp(self):
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.tol=0.001
    def tearDown(self):
        del self.domain
        del self.tol

class Test_DomainOnRipley(Test_Domain):
    def setUp(self):
        self.boundary_tag_list = [1, 2, 10, 20]
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain
        del self.boundary_tag_list

    def test_tagsContinuousFunction(self):
        ref_tags=[0]
        tags=ContinuousFunction(self.domain).getListOfTags()
        self.assertTrue(len(tags)==len(ref_tags), "tags list has wrong length.")
        for i in ref_tags: self.assertTrue(i in tags,"tag %s is missing."%i)
    def test_tagsFunction(self):
        ref_tags=[0]
        tags=Function(self.domain).getListOfTags()
        self.assertTrue(len(tags)==len(ref_tags), "tags list has wrong length.")
        for i in ref_tags: self.assertTrue(i in tags,"tag %s is missing."%i)
    def test_tagsReducedFunction(self):
        ref_tags=[0]
        tags=ReducedFunction(self.domain).getListOfTags()
        self.assertTrue(len(tags)==len(ref_tags), "tags list has wrong length.")
        for i in ref_tags: self.assertTrue(i in tags,"tag %s is missing."%i)
    def test_tagsFunctionOnBoundary(self):
        ref_tags=[1, 2, 10, 20]
        tags=FunctionOnBoundary(self.domain).getListOfTags()
        # For an MPI-distributed domain some tags may be missing
        if getMPISizeWorld() == 1: self.assertTrue(len(tags)==len(ref_tags), "tags list has wrong length.")
        for i in tags: self.assertTrue(i in ref_tags,"tag %s is missing."%i)
    def test_tagsReducedFunctionOnBoundary(self):
        ref_tags=[1, 2, 10, 20]
        tags=ReducedFunctionOnBoundary(self.domain).getListOfTags()
        # For an MPI-distributed domain some tags may be missing
        if getMPISizeWorld() == 1: self.assertTrue(len(tags)==len(ref_tags), "tags list has wrong length.")
        for i in tags: self.assertTrue(i in ref_tags,"tag %s is missing."%i)

class Test_DataOpsOnRipley(Test_Dump, Test_SetDataPointValue, Test_GlobalMinMax, Test_Lazy):
    def setUp(self):
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.domain_with_different_number_of_samples=Rectangle(n0=7*NE*NX-1, n1=3*NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.domain_with_different_number_of_data_points_per_sample=Rectangle(n0=7*NE*NX-1, n1=3*NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.domain_with_different_sample_ordering=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.filename_base=RIPLEY_WORKDIR
        self.mainfs=Function(self.domain)
        self.otherfs=Solution(self.domain)

    def tearDown(self):
        del self.domain
        del self.domain_with_different_number_of_samples
        del self.domain_with_different_number_of_data_points_per_sample
        del self.domain_with_different_sample_ordering
        del self.mainfs
        del self.otherfs

class Test_TableInterpolationOnRipley(Test_TableInterpolation):
    def setUp(self):
        self.domain = Brick(n0=NE*NXb-1, n1=NE*NYb-1, n2=NE*NZb-1, l0=1., l1=1., l2=1., d0=NXb, d1=NYb, d2=NZb)
        self.functionspaces=[ContinuousFunction(self.domain), Function(self.domain), ReducedFunction(self.domain),
            FunctionOnBoundary(self.domain), ReducedFunctionOnBoundary(self.domain)]
        #We aren't testing DiracDeltaFunctions
        self.xn=5 # number of grids on x axis
        self.yn=5 # number of grids on y axis
        self.zn=5

    def tearDown(self):
        del self.domain
        del self.functionspaces

class Test_CSVOnRipley(Test_saveCSV):
    def setUp(self):
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.linecount1=(NE*NX)*(NE*NY)+1  #see test_save1 for the meaning of these params
        self.linecount2=(NE*NX-1)*(NE*NY)+1
        self.line_expected=[1./(NE*NX-1), 0., 1./(NE*NX-1)]

    def tearDown(self):
        del self.domain

    #This test checks to see that all FunctionSpaces can be saved
    def test_singleFS(self):
        fname=os.path.join(RIPLEY_WORKDIR, "test_singlefs.csv")
        fss=[ContinuousFunction(self.domain), Function(self.domain), ReducedFunction(self.domain),
        FunctionOnBoundary(self.domain), ReducedFunctionOnBoundary(self.domain),
        DiracDeltaFunctions(self.domain)]
        for f in fss:
            d=Data(7,f)
            print("Testing "+str(f)+"\n")
            saveDataCSV(fname, D=d)

    def test_multiFS(self):
        fname=os.path.join(RIPLEY_WORKDIR, "test_multifs.csv")
        sol=Data(8,Solution(self.domain))
        ctsfn=Data(9,ContinuousFunction(self.domain))
        #test line 0
        dirac=Data(-1,DiracDeltaFunctions(self.domain))
        saveDataCSV(fname, A=sol, B=ctsfn, C=dirac)
        #test line 1
        fun=Data(5,Function(self.domain))
        rfun=Data(3,ReducedFunction(self.domain))
        saveDataCSV(fname, A=sol,B=ctsfn,C=fun, D=rfun)
        #test line 2
        bound=Data(1,FunctionOnBoundary(self.domain))
        rbound=Data(3,ReducedFunctionOnBoundary(self.domain))
        saveDataCSV(fname,A=sol,B=ctsfn,C=bound, D=rbound)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_SharedOnRipley))
    suite.addTest(unittest.makeSuite(Test_DataOpsOnRipley))
    suite.addTest(unittest.makeSuite(Test_DomainOnRipley))
    suite.addTest(unittest.makeSuite(Test_TableInterpolationOnRipley))
    suite.addTest(unittest.makeSuite(Test_CSVOnRipley))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

