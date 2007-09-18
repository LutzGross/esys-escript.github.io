#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import unittest
import tempfile

from esys.escript import *
from esys.finley import Rectangle
import sys
import os

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=FINLEY_TEST_DATA+"/data_meshes/"


NE=4 # number elements, must be even
class Test_Dump(unittest.TestCase):
   arg0=9.81
   arg1=numarray.array([3.098, -3.111])
   arg2=numarray.array([[3.82, -3.81, -0.957, 0.892, -1.367], [-4.589, -1.835, -2.679, -1.517, -4.2515], [-4.909, 1.634, -2.883, -2.135, 1.187], [0.6431, 4.638, -4.616, -0.196, -4.370]])
   arg3=numarray.array([[[-2.3667, -0.040], [-4.7398, -3.2412]], [[-2.125, -2.240], [2.237, -4.279]], [[0.68720, 2.4059], [-2.4964, 3.17453]], [[-4.907, -4.9431], [-0.3604, 0.4269]], [[1.4179, 3.326], [1.356, -0.4610]], [[3.378, 2.0902], [-2.6857, 1.3585]]])
   arg4=numarray.array([[[[-3.810, -1.3597, -1.5307, 1.099], [-1.828, 0.2526, -1.4429, 2.326], [4.9732, -2.063, 1.3153, -3.809]], [[-4.8902, -4.714, 1.520, -1.931], [-3.8847, 4.3867, 1.894030, 2.432], [-1.2082, -0.8304, 2.2612, 4.6399]]], [[[-4.5922, -3.309, -0.8171, -0.7210], [2.8051, -4.93047, 0.08450, 4.3824], [0.43204, 2.1908, 4.512633, -1.8218]], [[2.2493, -4.190, -2.3893, -4.147], [-2.104, -4.635, -4.2767, -3.53151], [-2.351, -1.6614, 2.9385, 4.099]]], [[[1.710, 0.2235, -3.4917, 0.8713], [-0.2881, 4.6278, 3.603, -2.1211], [-0.565, 4.294, -2.210827, -0.37651]], [[0.6578, -2.869, -2.490, -4.789], [3.232, 2.483, 0.9531, 2.260], [-1.785, 0.42156, -1.8379, 4.212]]]])

   def _diffDataObjects(self,d_ref,file):
       d_ref.dump(file)
       d=load(file, self.domain)
       self.failUnless(not d.isEmpty(),"data in %s are empty."%file)
       self.failUnless(d_ref.getFunctionSpace() == d.getFunctionSpace(), "wrong function space in %s."%file)
       self.failUnless(d_ref.getRank() == d.getRank(), "different rank in %s. "%file)
       self.failUnless(d_ref.getShape() == d.getShape(), "different shape %s. "%file)
       self.failUnless(Lsup(d_ref-d)<=0., "different entries %s."%file)

   #===========================================================================
   def test_DumpAndLoad_Constant_Solution_Rank0(self):
       file="constant_solution_rank0.nc"
       d=Data(self.arg0,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank1(self):
       file="constant_solution_rank1.nc"
       d=Data(self.arg1,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank2(self):
       file="constant_solution_rank2.nc"
       d=Data(self.arg2,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank3(self):
       file="constant_solution_rank3.nc"
       d=Data(self.arg3,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank4(self):
       file="constant_solution_rank4.nc"
       d=Data(self.arg4,self.solution)
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedSolution_Rank0(self):
       file="constant_reduced_solution_rank0.nc"
       d=Data(self.arg0,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank1(self):
       file="constant_reduced_solution_rank1.nc"
       d=Data(self.arg1,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank2(self):
       file="constant_reduced_solution_rank2.nc"
       d=Data(self.arg2,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank3(self):
       file="constant_reduced_solution_rank3.nc"
       d=Data(self.arg3,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank4(self):
       file="constant_reduced_solution_rank4.nc"
       d=Data(self.arg4,self.reduced_solution)
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Constant_ContinuousFunction_Rank0(self):
       file="constant_continuous_function_rank0.nc"
       d=Data(self.arg0,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank1(self):
       file="constant_continuous_function_rank1.nc"
       d=Data(self.arg1,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank2(self):
       file="constant_continuous_function_rank2.nc"
       d=Data(self.arg2,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank3(self):
       file="constant_continuous_function_rank3.nc"
       d=Data(self.arg3,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank4(self):
       file="constant_continuous_function_rank4.nc"
       d=Data(self.arg4,self.continuous_function)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Constant_Function_Rank0(self):
       file="constant_function_rank0.nc"
       d=Data(self.arg0,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank1(self):
       file="constant_function_rank1.nc"
       d=Data(self.arg1,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank2(self):
       file="constant_function_rank2.nc"
       d=Data(self.arg2,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank3(self):
       file="constant_function_rank3.nc"
       d=Data(self.arg3,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank4(self):
       file="constant_function_rank4.nc"
       d=Data(self.arg4,self.function)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank0(self):
       file="constant_function_on_boundary_rank0.nc"
       d=Data(self.arg0,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank1(self):
       file="constant_function_on_boundary_rank1.nc"
       d=Data(self.arg1,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank2(self):
       file="constant_function_on_boundary_rank2.nc"
       d=Data(self.arg2,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank3(self):
       file="constant_function_on_boundary_rank3.nc"
       d=Data(self.arg3,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank4(self):
       file="constant_function_on_boundary_rank4.nc"
       d=Data(self.arg4,self.function_on_boundary)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Solution_Rank0(self):
       file="expanded_solution_rank0.nc"
       d=Data(length(self.solution.getX())*self.arg0,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank1(self):
       file="expanded_solution_rank1.nc"
       d=Data(length(self.solution.getX())*self.arg1,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank2(self):
       file="expanded_solution_rank2.nc"
       d=Data(length(self.solution.getX())*self.arg2,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank3(self):
       file="expanded_solution_rank3.nc"
       d=Data(length(self.solution.getX())*self.arg3,self.solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank4(self):
       file="expanded_solution_rank4.nc"
       d=Data(length(self.solution.getX())*self.arg4,self.solution)
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedSolution_Rank0(self):
       file="expanded_reduced_solution_rank0.nc"
       d=Data(length(self.reduced_solution.getX())*self.arg0,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank1(self):
       file="expanded_reduced_solution_rank1.nc"
       d=Data(length(self.reduced_solution.getX())*self.arg1,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank2(self):
       file="expanded_reduced_solution_rank2.nc"
       d=Data(length(self.reduced_solution.getX())*self.arg2,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank3(self):
       file="expanded_reduced_solution_rank3.nc"
       d=Data(length(self.reduced_solution.getX())*self.arg3,self.reduced_solution)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank4(self):
       file="expanded_reduced_solution_rank4.nc"
       d=Data(length(self.reduced_solution.getX())*self.arg4,self.reduced_solution)
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank0(self):
       file="expanded_continuous_function_rank0.nc"
       d=Data(length(self.continuous_function.getX())*self.arg0,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank1(self):
       file="expanded_continuous_function_rank1.nc"
       d=Data(length(self.continuous_function.getX())*self.arg1,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank2(self):
       file="expanded_continuous_function_rank2.nc"
       d=Data(length(self.continuous_function.getX())*self.arg2,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank3(self):
       file="expanded_continuous_function_rank3.nc"
       d=Data(length(self.continuous_function.getX())*self.arg3,self.continuous_function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank4(self):
       file="expanded_continuous_function_rank4.nc"
       d=Data(length(self.continuous_function.getX())*self.arg4,self.continuous_function)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Function_Rank0(self):
       file="expanded_function_rank0.nc"
       d=Data(length(self.function.getX())*self.arg0,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank1(self):
       file="expanded_function_rank1.nc"
       d=Data(length(self.function.getX())*self.arg1,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank2(self):
       file="expanded_function_rank2.nc"
       d=Data(length(self.function.getX())*self.arg2,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank3(self):
       file="expanded_function_rank3.nc"
       d=Data(length(self.function.getX())*self.arg3,self.function)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank4(self):
       file="expanded_function_rank4.nc"
       d=Data(length(self.function.getX())*self.arg4,self.function)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank0(self):
       file="expanded_function_on_boundary_rank0.nc"
       d=Data(length(self.function_on_boundary.getX())*self.arg0,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank1(self):
       file="expanded_function_on_boundary_rank1.nc"
       d=Data(length(self.function_on_boundary.getX())*self.arg1,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank2(self):
       file="expanded_function_on_boundary_rank2.nc"
       d=Data(length(self.function_on_boundary.getX())*self.arg2,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank3(self):
       file="expanded_function_on_boundary_rank3.nc"
       d=Data(length(self.function_on_boundary.getX())*self.arg3,self.function_on_boundary)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank4(self):
       file="expanded_function_on_boundary_rank4.nc"
       d=Data(length(self.function_on_boundary.getX())*self.arg4,self.function_on_boundary)
       self._diffDataObjects(d,file)



class Test_DumpOnFinley(Test_Dump):
   def setUp(self):
       self.domain =Rectangle(NE,NE+1,2)
       self.solution = Solution(self.domain) 
       self.reduced_solution = ReducedSolution(self.domain) 
       self.continuous_function = ContinuousFunction(self.domain) 
       self.function = Function(self.domain) 
       self.function_on_boundary = FunctionOnBoundary(self.domain) 

   def tearDown(self):
       del self.function
       del self.function_on_boundary
       del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_DumpOnFinley))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
   
