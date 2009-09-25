
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
import tempfile

from esys.escript import *
from esys.finley import Rectangle
import sys
import os
from test_objects import Test_Dump, Test_SetDataPointValue, Test_saveCSV, Test_TableInterpolation
from test_objects import Test_Domain, Test_GlobalMinMax

from test_shared import Test_Shared

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

NE=4 # number elements, must be even

class Test_SharedOnFinley(Test_Shared):
  def setUp(self):
	self.domain=Rectangle(NE,NE)
	self.tol=0.001
  def tearDown(self):
	del self.domain
	del self.tol

class Test_DomainOnFinley(Test_Domain):
   def setUp(self):
       self.boundary_tag_list = [1, 2, 10, 20]
       self.domain =Rectangle(NE,NE+1,2)
   def tearDown(self):
       del self.domain
       del self.boundary_tag_list 

   def test_tagsContinuousFunction(self):
       ref_tags=[0]
       tags=ContinuousFunction(self.domain).getListOfTags()
       self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in ref_tags: self.failUnless(i in tags,"tag %s is missing."%i)

   def test_tagsFunction(self):
       ref_tags=[0]
       tags=Function(self.domain).getListOfTags()
       self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in ref_tags: self.failUnless(i in tags,"tag %s is missing."%i)
   def test_tagsReducedFunction(self):
       ref_tags=[0]
       tags=ReducedFunction(self.domain).getListOfTags()
       self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in ref_tags: self.failUnless(i in tags,"tag %s is missing."%i)
   def test_tagsFunctionOnBoundary(self):
       ref_tags=[1, 2, 10, 20]
       tags=FunctionOnBoundary(self.domain).getListOfTags()
       # For an MPI-distributed domain some tags may be missing
       if getMPISizeWorld() == 1: self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in tags: self.failUnless(i in ref_tags,"tag %s is missing."%i)
   def test_tagsReducedFunctionOnBoundary(self):
       ref_tags=[1, 2, 10, 20]
       tags=ReducedFunctionOnBoundary(self.domain).getListOfTags()
       # For an MPI-distributed domain some tags may be missing
       if getMPISizeWorld() == 1: self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in tags: self.failUnless(i in ref_tags,"tag %s is missing."%i)
   def test_tagsFunctionOnContactOne(self):
       ref_tags=[]
       tags=FunctionOnContactOne(self.domain).getListOfTags()
       self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in ref_tags: self.failUnless(i in tags,"tag %s is missing."%i)
   def test_tagsFunctionOnContactZero(self):
       ref_tags=[]
       tags=FunctionOnContactZero(self.domain).getListOfTags()
       self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in ref_tags: self.failUnless(i in tags,"tag %s is missing."%i)
   def test_tagsReducedFunctionOnContactOne(self):
       ref_tags=[]
       tags=ReducedFunctionOnContactOne(self.domain).getListOfTags()
       self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in ref_tags: self.failUnless(i in tags,"tag %s is missing."%i)
   def test_tagsReducedFunctionOnContactZero(self):
       ref_tags=[]
       tags=ReducedFunctionOnContactZero(self.domain).getListOfTags()
       self.failUnless(len(tags)==len(ref_tags), "tags list has wrong length.")
       for i in ref_tags: self.failUnless(i in tags,"tag %s is missing."%i)

class Test_DataOpsOnFinley(Test_Dump, Test_SetDataPointValue, Test_GlobalMinMax):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
       self.domain_with_different_number_of_samples =Rectangle(2*NE,NE+1,2)
       self.domain_with_different_number_of_data_points_per_sample =Rectangle(2*NE,NE+1,2,integrationOrder=2)
       self.domain_with_different_sample_ordering =Rectangle(NE,NE+1,2, optimize=True)
       self.filename_base=FINLEY_WORKDIR

   def tearDown(self):
       del self.domain
       del self.domain_with_different_number_of_samples
       del self.domain_with_different_number_of_data_points_per_sample
       del self.domain_with_different_sample_ordering
       


class Test_TableInterpolationOnFinley(Test_TableInterpolation):
    def setUp(self):
	self.domain=Rectangle(4,4)
	self.functionspaces=[ContinuousFunction(self.domain), Function(self.domain), ReducedFunction(self.domain),
	    FunctionOnBoundary(self.domain), ReducedFunctionOnBoundary(self.domain), 
	    FunctionOnContactZero(self.domain), FunctionOnContactOne(self.domain),
	    ReducedFunctionOnContactZero(self.domain), ReducedFunctionOnContactOne(self.domain)]
	    #We aren't testing DiracDeltaFunction
	self.xn=3	# number of grids on x axis
	self.yn=3	# number of grids on y axis

    def tearDown(self):
	del self.domain
	del self.functionspaces

	
		
		
class Test_CSVOnFinley(Test_saveCSV):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
       self.linecount1=20		#see test_save1 for the meaning of these params
       self.linecount2=69
       
   def tearDown(self):
       del self.domain
       
   #This test checks to see that all FunctionSpaces can be saved
   def test_singleFS(self):
	fname="test_singlefs.csv"
	fss=[ContinuousFunction(self.domain), Function(self.domain), ReducedFunction(self.domain),
	FunctionOnBoundary(self.domain), ReducedFunctionOnBoundary(self.domain), 
	FunctionOnContactZero(self.domain), FunctionOnContactOne(self.domain),
	ReducedFunctionOnContactZero(self.domain), ReducedFunctionOnContactOne(self.domain),
	DiracDeltaFunction(self.domain)]
	for f in fss:
		d=Data(7,f)
		print "Testing "+str(f)+"\n"
		saveDataCSV(fname, D=d)

   def test_multiFS(self):
	fname="test_multifs.csv"
	sol=Data(8,Solution(self.domain))
	ctsfn=Data(9,ContinuousFunction(self.domain))
	#test line 0
	dirac=Data(-1,DiracDeltaFunction(self.domain))
	saveDataCSV(fname, A=sol, B=ctsfn, C=dirac)
	#test line 1
	fun=Data(5,Function(self.domain))
	rfun=Data(3,ReducedFunction(self.domain))
	saveDataCSV(fname, A=sol,B=ctsfn,C=fun, D=rfun)
	#test line 2
	bound=Data(1,FunctionOnBoundary(self.domain))
	rbound=Data(3,ReducedFunctionOnBoundary(self.domain))
	saveDataCSV(fname,A=sol,B=ctsfn,C=bound, D=rbound)
	#test line 3
	conzz=Data(7,FunctionOnContactZero(self.domain))
	rconz=Data(8,ReducedFunctionOnContactZero(self.domain))
	saveDataCSV(fname,A=sol,B=ctsfn, C=conzz, D=rconz)
	#check for cross line exceptions
	self.failUnlessRaises(RuntimeError, saveDataCSV, fname, A=dirac, B=rfun)
	self.failUnlessRaises(RuntimeError, saveDataCSV, fname, A=bound, B=conzz)

	
if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_SharedOnFinley))
   suite.addTest(unittest.makeSuite(Test_DataOpsOnFinley))
   suite.addTest(unittest.makeSuite(Test_DomainOnFinley))
   suite.addTest(unittest.makeSuite(Test_TableInterpolationOnFinley))
   suite.addTest(unittest.makeSuite(Test_CSVOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

