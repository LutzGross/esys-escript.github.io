
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

import os
import sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.dudley import Rectangle, Brick
from test_objects import Test_Dump, Test_SetDataPointValue, Test_saveCSV, \
        Test_TableInterpolation, Test_Domain, Test_Lazy
from test_shared import Test_Shared

try:
     DUDLEY_WORKDIR=os.environ['DUDLEY_WORKDIR']
except KeyError:
     DUDLEY_WORKDIR='.'

NE = max(4, getMPISizeWorld())

class Test_SharedOnDudley(Test_Shared):
  def setUp(self):
        self.domain=Rectangle(NE,NE)
        self.tol=0.001
  def tearDown(self):
        del self.domain
        del self.tol

class Test_DomainOnDudley(Test_Domain):
   def setUp(self):
       self.boundary_tag_list = [1, 2, 10, 20]
       self.domain = Rectangle(NE,NE+1)
       self.rdomain=self.domain

   def tearDown(self):
       del self.domain
       del self.rdomain
       del self.boundary_tag_list

   def test_setXError(self):
       domain=Rectangle(NE,NE)
       x=domain.getX()
       z=interpolate(x, Function(domain))
       self.assertRaises(ValueError, domain.setX, z)
       del x
       del z
       del domain

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

class Test_DumpOnDudley(Test_Dump):
   def setUp(self):
       self.domain = Rectangle(NE,NE+1)
       self.domain_with_different_number_of_samples = Rectangle(2*NE,NE+1)
       self.domain_with_different_number_of_data_points_per_sample = Rectangle(2*NE,NE+1,integrationOrder=2)
       self.domain_with_different_sample_ordering = Rectangle(NE, NE+1, optimize=True)
       self.filename_base = DUDLEY_WORKDIR

   def tearDown(self):
       del self.domain
       del self.domain_with_different_number_of_samples
       del self.domain_with_different_number_of_data_points_per_sample
       del self.domain_with_different_sample_ordering

class Test_SetDataPointValueOnDudley(Test_SetDataPointValue):
   def setUp(self):
       self.domain = Rectangle(NE,NE+1)

   def tearDown(self):
       del self.domain

class Test_LazyOnDudley(Test_Lazy):
   def setUp(self):
       self.domain = Rectangle(NE,NE+1)
       self.mainfs = Function(self.domain)
       self.otherfs = Solution(self.domain)

   def tearDown(self):
       del self.domain
       del self.mainfs
       del self.otherfs

class Test_TableInterpolationOnDudley(Test_TableInterpolation):
    def setUp(self):
        self.domain=Brick(4,4,4)
        self.functionspaces=[ContinuousFunction(self.domain), Function(self.domain), ReducedFunction(self.domain),
            FunctionOnBoundary(self.domain), ReducedFunctionOnBoundary(self.domain)]
            #We aren't testing DiracDeltaFunctions
        self.xn=5       # number of grids on x axis
        self.yn=5       # number of grids on y axis
        self.zn=5

    def tearDown(self):
        del self.domain
        del self.functionspaces


class Test_CSVOnDudley(Test_saveCSV):
    def setUp(self):
        self.workdir=DUDLEY_WORKDIR
        NE0=NE
        NE1=NE+1
        self.domain=Rectangle(NE0,NE1)
        self.functionspaces=[ ContinuousFunction ]
        # number of total data points for each function space
        self.linecounts=[ (NE0+1)*(NE1+1)+1 ]
        # number of masked points, i.e. where X[0] is non-zero
        self.linecounts_masked=[ NE0*(NE1+1)+1 ]
        # expected values in first line of masked data = [ X[:], X[0] ]
        self.firstline=[ [1./NE0, 0., 1./NE0] ]

        if getMPISizeWorld() == 1:
            self.functionspaces += [ Function, ReducedFunction,
                    FunctionOnBoundary, ReducedFunctionOnBoundary ]
            self.linecounts += [ 121, 41, 37, 19 ]
            self.linecounts_masked += [ 116, 41, 27, 14 ]
            self.firstline += [ 
                    [.125, 0., .125],
                    [.16666666666666667, .0666666666666667, .166666666666667],
                    [.05283121635129676, 0., .05283121635129676],
                    [.125, 0., .125]
                    ]

    def tearDown(self):
        del self.domain

    @unittest.skipIf(getMPISizeWorld() > 1, "Skipping since MPI size > 1")
    def test_csv_multiFS(self):
        fname=os.path.join(DUDLEY_WORKDIR, "test_multifs.csv")
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
    run_tests(__name__, exit_on_failure=True)

