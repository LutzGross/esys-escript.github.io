# $Id:$
"""
Test suite for data objects. at the moment for dump and load only.

The tests must be linked with some function space class object in the setUp method:
to run the use:

   from esys.finley import Brick
   class Test_DumpOnFinley(Test_Dump):
       def setUp(self):
          self.domain =Rectangle(NE,NE+1,2)
          self.domain_with_different_number_of_samples =Rectangle(2*NE,NE+1,2)
          self.domain_with_different_number_of_data_points_per_sample =Rectangle(2*NE,NE+1,2,integrationOrder=2)
          self.domain_with_different_sample_ordering =Rectangle(1,(NE+1)*NE,2)
          self.filebase="."

   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_DumpOnFinley))
   unittest.TextTestRunner(verbosity=2).run(suite)

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2006, 2007 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision:$"
__date__="$Date:$"

import unittest
import os
import numarray
from esys.escript import *

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
       file=os.path.join(self.filebase,"constant_solution_rank0.nc")
       d=Data(self.arg0,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank1(self):
       file=os.path.join(self.filebase,"constant_solution_rank1.nc")
       d=Data(self.arg1,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank2(self):
       file=os.path.join(self.filebase,"constant_solution_rank2.nc")
       d=Data(self.arg2,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank3(self):
       file=os.path.join(self.filebase,"constant_solution_rank3.nc")
       d=Data(self.arg3,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Solution_Rank4(self):
       file=os.path.join(self.filebase,"constant_solution_rank4.nc")
       d=Data(self.arg4,Solution(self.domain))
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedSolution_Rank0(self):
       file=os.path.join(self.filebase,"constant_reduced_solution_rank0.nc")
       d=Data(self.arg0,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank1(self):
       file=os.path.join(self.filebase,"constant_reduced_solution_rank1.nc")
       d=Data(self.arg1,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank2(self):
       file=os.path.join(self.filebase,"constant_reduced_solution_rank2.nc")
       d=Data(self.arg2,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank3(self):
       file=os.path.join(self.filebase,"constant_reduced_solution_rank3.nc")
       d=Data(self.arg3,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank4(self):
       file=os.path.join(self.filebase,"constant_reduced_solution_rank4.nc")
       d=Data(self.arg4,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Constant_ContinuousFunction_Rank0(self):
       file=os.path.join(self.filebase,"constant_continuous_function_rank0.nc")
       d=Data(self.arg0,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank1(self):
       file=os.path.join(self.filebase,"constant_continuous_function_rank1.nc")
       d=Data(self.arg1,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank2(self):
       file=os.path.join(self.filebase,"constant_continuous_function_rank2.nc")
       d=Data(self.arg2,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank3(self):
       file=os.path.join(self.filebase,"constant_continuous_function_rank3.nc")
       d=Data(self.arg3,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank4(self):
       file=os.path.join(self.filebase,"constant_continuous_function_rank4.nc")
       d=Data(self.arg4,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Constant_Function_Rank0(self):
       file=os.path.join(self.filebase,"constant_function_rank0.nc")
       d=Data(self.arg0,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank1(self):
       file=os.path.join(self.filebase,"constant_function_rank1.nc")
       d=Data(self.arg1,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank2(self):
       file=os.path.join(self.filebase,"constant_function_rank2.nc")
       d=Data(self.arg2,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank3(self):
       file=os.path.join(self.filebase,"constant_function_rank3.nc")
       d=Data(self.arg3,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_Function_Rank4(self):
       file=os.path.join(self.filebase,"constant_function_rank4.nc")
       d=Data(self.arg4,Function(self.domain))
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank0(self):
       file=os.path.join(self.filebase,"constant_function_on_boundary_rank0.nc")
       d=Data(self.arg0,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank1(self):
       file=os.path.join(self.filebase,"constant_function_on_boundary_rank1.nc")
       d=Data(self.arg1,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank2(self):
       file=os.path.join(self.filebase,"constant_function_on_boundary_rank2.nc")
       d=Data(self.arg2,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank3(self):
       file=os.path.join(self.filebase,"constant_function_on_boundary_rank3.nc")
       d=Data(self.arg3,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank4(self):
       file=os.path.join(self.filebase,"constant_function_on_boundary_rank4.nc")
       d=Data(self.arg4,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Solution_Rank0(self):
       file=os.path.join(self.filebase,"expanded_solution_rank0.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Solution(self.domain).getX())*self.arg0,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank1(self):
       file=os.path.join(self.filebase,"expanded_solution_rank1.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Solution(self.domain).getX())*self.arg1,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank2(self):
       file=os.path.join(self.filebase,"expanded_solution_rank2.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Solution(self.domain).getX())*self.arg2,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank3(self):
       file=os.path.join(self.filebase,"expanded_solution_rank3.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Solution(self.domain).getX())*self.arg3,Solution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Solution_Rank4(self):
       file=os.path.join(self.filebase,"expanded_solution_rank4.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Solution(self.domain).getX())*self.arg4,Solution(self.domain))
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedSolution_Rank0(self):
       file=os.path.join(self.filebase,"expanded_reduced_solution_rank0.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ReducedSolution(self.domain).getX())*self.arg0,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank1(self):
       file=os.path.join(self.filebase,"expanded_reduced_solution_rank1.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ReducedSolution(self.domain).getX())*self.arg1,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank2(self):
       file=os.path.join(self.filebase,"expanded_reduced_solution_rank2.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ReducedSolution(self.domain).getX())*self.arg2,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank3(self):
       file=os.path.join(self.filebase,"expanded_reduced_solution_rank3.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ReducedSolution(self.domain).getX())*self.arg3,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank4(self):
       file=os.path.join(self.filebase,"expanded_reduced_solution_rank4.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ReducedSolution(self.domain).getX())*self.arg4,ReducedSolution(self.domain))
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank0(self):
       file=os.path.join(self.filebase,"expanded_continuous_function_rank0.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ContinuousFunction(self.domain).getX())*self.arg0,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank1(self):
       file=os.path.join(self.filebase,"expanded_continuous_function_rank1.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ContinuousFunction(self.domain).getX())*self.arg1,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank2(self):
       file=os.path.join(self.filebase,"expanded_continuous_function_rank2.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ContinuousFunction(self.domain).getX())*self.arg2,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank3(self):
       file=os.path.join(self.filebase,"expanded_continuous_function_rank3.nc")
       d=Data(length(ContinuousFunction(self.domain).getX())*self.arg3,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank4(self):
       file=os.path.join(self.filebase,"expanded_continuous_function_rank4.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(ContinuousFunction(self.domain).getX())*self.arg4,ContinuousFunction(self.domain))
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Function_Rank0(self):
       file=os.path.join(self.filebase,"expanded_function_rank0.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       # elements are not in different order: self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Function(self.domain).getX())*self.arg0,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank1(self):
       file=os.path.join(self.filebase,"expanded_function_rank1.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       # elements are not in different order: self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Function(self.domain).getX())*self.arg1,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank2(self):
       file=os.path.join(self.filebase,"expanded_function_rank2.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       # elements are not in different order: self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Function(self.domain).getX())*self.arg2,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank3(self):
       file=os.path.join(self.filebase,"expanded_function_rank3.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       # elements are not in different order: self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Function(self.domain).getX())*self.arg3,Function(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_Function_Rank4(self):
       file=os.path.join(self.filebase,"expanded_function_rank4.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       # elements are not in different order: self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(Function(self.domain).getX())*self.arg4,Function(self.domain))
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank0(self):
       file=os.path.join(self.filebase,"expanded_function_on_boundary_rank0.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg0,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank1(self):
       file=os.path.join(self.filebase,"expanded_function_on_boundary_rank1.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg1,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank2(self):
       file=os.path.join(self.filebase,"expanded_function_on_boundary_rank2.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg2,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank3(self):
       file=os.path.join(self.filebase,"expanded_function_on_boundary_rank3.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg3,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank4(self):
       file=os.path.join(self.filebase,"expanded_function_on_boundary_rank4.nc")
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_samples)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_number_of_data_points_per_sample)
       self.failUnlessRaises(RuntimeError, load, file, self.domain_with_different_sample_ordering)
       d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg4,FunctionOnBoundary(self.domain))
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Tagged_Solution_Rank0(self):
       file=os.path.join(self.filebase,"tagged_solution_rank0.nc")
       d=Data(self.arg0,Solution(self.domain))
       d.setTaggedValue(1,self.arg0*2)
       d.setTaggedValue(10,self.arg0*3)
       d.setTaggedValue(100,self.arg0*4)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Tagged_Solution_Rank0(self):
       file=os.path.join(self.filebase,"tagged_solution_rank0.nc")
       d=Data(self.arg0,Solution(self.domain))
       d.setTaggedValue(1,self.arg0*2)
       d.setTaggedValue(10,self.arg0*3)
       d.setTaggedValue(100,self.arg0*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Solution_Rank1(self):
       file=os.path.join(self.filebase,"tagged_solution_rank1.nc")
       d=Data(self.arg1,Solution(self.domain))
       d.setTaggedValue(1,self.arg1*2)
       d.setTaggedValue(10,self.arg1*3)
       d.setTaggedValue(100,self.arg1*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Solution_Rank2(self):
       file=os.path.join(self.filebase,"tagged_solution_rank2.nc")
       d=Data(self.arg2,Solution(self.domain))
       d.setTaggedValue(1,self.arg2*2)
       d.setTaggedValue(10,self.arg2*3)
       d.setTaggedValue(100,self.arg2*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Solution_Rank3(self):
       file=os.path.join(self.filebase,"tagged_solution_rank3.nc")
       d=Data(self.arg3,Solution(self.domain))
       d.setTaggedValue(1,self.arg3*2)
       d.setTaggedValue(10,self.arg3*3)
       d.setTaggedValue(100,self.arg3*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Solution_Rank4(self):
       file=os.path.join(self.filebase,"tagged_solution_rank4.nc")
       d=Data(self.arg4,Solution(self.domain))
       d.setTaggedValue(1,self.arg4*2)
       d.setTaggedValue(10,self.arg4*3)
       d.setTaggedValue(100,self.arg4*4)
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Tagged_ReducedSolution_Rank0(self):
       file=os.path.join(self.filebase,"tagged_reduced_solution_rank0.nc")
       d=Data(self.arg0,ReducedSolution(self.domain))
       d.setTaggedValue(1,self.arg0*2)
       d.setTaggedValue(10,self.arg0*3)
       d.setTaggedValue(100,self.arg0*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ReducedSolution_Rank1(self):
       file=os.path.join(self.filebase,"tagged_reduced_solution_rank1.nc")
       d=Data(self.arg1,ReducedSolution(self.domain))
       d.setTaggedValue(1,self.arg1*2)
       d.setTaggedValue(10,self.arg1*3)
       d.setTaggedValue(100,self.arg1*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ReducedSolution_Rank2(self):
       file=os.path.join(self.filebase,"tagged_reduced_solution_rank2.nc")
       d=Data(self.arg2,ReducedSolution(self.domain))
       d.setTaggedValue(1,self.arg2*2)
       d.setTaggedValue(10,self.arg2*3)
       d.setTaggedValue(100,self.arg2*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ReducedSolution_Rank3(self):
       file=os.path.join(self.filebase,"tagged_reduced_solution_rank3.nc")
       d=Data(self.arg3,ReducedSolution(self.domain))
       d.setTaggedValue(1,self.arg3*2)
       d.setTaggedValue(10,self.arg3*3)
       d.setTaggedValue(100,self.arg3*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ReducedSolution_Rank4(self):
       file=os.path.join(self.filebase,"tagged_reduced_solution_rank4.nc")
       d=Data(self.arg4,ReducedSolution(self.domain))
       d.setTaggedValue(1,self.arg4*2)
       d.setTaggedValue(10,self.arg4*3)
       d.setTaggedValue(100,self.arg4*4)
       self._diffDataObjects(d,file)
   #===========================================================================
   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank0(self):
       file=os.path.join(self.filebase,"tagged_continuous_function_rank0.nc")
       d=Data(self.arg0,ContinuousFunction(self.domain))
       d.setTaggedValue(1,self.arg0*2)
       d.setTaggedValue(10,self.arg0*3)
       d.setTaggedValue(100,self.arg0*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank1(self):
       file=os.path.join(self.filebase,"tagged_continuous_function_rank1.nc")
       d=Data(self.arg1,ContinuousFunction(self.domain))
       d.setTaggedValue(1,self.arg1*2)
       d.setTaggedValue(10,self.arg1*3)
       d.setTaggedValue(100,self.arg1*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank2(self):
       file=os.path.join(self.filebase,"tagged_continuous_function_rank2.nc")
       d=Data(self.arg2,ContinuousFunction(self.domain))
       d.setTaggedValue(1,self.arg2*2)
       d.setTaggedValue(10,self.arg2*3)
       d.setTaggedValue(100,self.arg2*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank3(self):
       file=os.path.join(self.filebase,"tagged_continuous_function_rank3.nc")
       d=Data(self.arg3,ContinuousFunction(self.domain))
       d.setTaggedValue(1,self.arg3*2)
       d.setTaggedValue(10,self.arg3*3)
       d.setTaggedValue(100,self.arg3*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank4(self):
       file=os.path.join(self.filebase,"tagged_continuous_function_rank4.nc")
       d=Data(self.arg4,ContinuousFunction(self.domain))
       d.setTaggedValue(1,self.arg4*2)
       d.setTaggedValue(10,self.arg4*3)
       d.setTaggedValue(100,self.arg4*4)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Tagged_Function_Rank0(self):
       file=os.path.join(self.filebase,"tagged_function_rank0.nc")
       d=Data(self.arg0,Function(self.domain))
       d.setTaggedValue(1,self.arg0*2)
       d.setTaggedValue(10,self.arg0*3)
       d.setTaggedValue(100,self.arg0*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Function_Rank1(self):
       file=os.path.join(self.filebase,"tagged_function_rank1.nc")
       d=Data(self.arg1,Function(self.domain))
       d.setTaggedValue(1,self.arg1*2)
       d.setTaggedValue(10,self.arg1*3)
       d.setTaggedValue(100,self.arg1*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Function_Rank2(self):
       file=os.path.join(self.filebase,"tagged_function_rank2.nc")
       d=Data(self.arg2,Function(self.domain))
       d.setTaggedValue(1,self.arg2*2)
       d.setTaggedValue(10,self.arg2*3)
       d.setTaggedValue(100,self.arg2*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Function_Rank3(self):
       file=os.path.join(self.filebase,"tagged_function_rank3.nc")
       d=Data(self.arg3,Function(self.domain))
       d.setTaggedValue(1,self.arg3*2)
       d.setTaggedValue(10,self.arg3*3)
       d.setTaggedValue(100,self.arg3*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_Function_Rank4(self):
       file=os.path.join(self.filebase,"tagged_function_rank4.nc")
       d=Data(self.arg4,Function(self.domain))
       d.setTaggedValue(1,self.arg4*2)
       d.setTaggedValue(10,self.arg4*3)
       d.setTaggedValue(100,self.arg4*4)
       self._diffDataObjects(d,file)

   #===========================================================================
   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank0(self):
       file=os.path.join(self.filebase,"tagged_function_on_boundary_rank0.nc")
       d=Data(self.arg0,FunctionOnBoundary(self.domain))
       d.setTaggedValue(1,self.arg0*2)
       d.setTaggedValue(10,self.arg0*3)
       d.setTaggedValue(100,self.arg0*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank1(self):
       file=os.path.join(self.filebase,"tagged_function_on_boundary_rank1.nc")
       d=Data(self.arg1,FunctionOnBoundary(self.domain))
       d.setTaggedValue(1,self.arg1*2)
       d.setTaggedValue(10,self.arg1*3)
       d.setTaggedValue(100,self.arg1*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank2(self):
       file=os.path.join(self.filebase,"tagged_function_on_boundary_rank2.nc")
       d=Data(self.arg2,FunctionOnBoundary(self.domain))
       d.setTaggedValue(1,self.arg2*2)
       d.setTaggedValue(10,self.arg2*3)
       d.setTaggedValue(100,self.arg2*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank3(self):
       file=os.path.join(self.filebase,"tagged_function_on_boundary_rank3.nc")
       d=Data(self.arg3,FunctionOnBoundary(self.domain))
       d.setTaggedValue(1,self.arg3*2)
       d.setTaggedValue(10,self.arg3*3)
       d.setTaggedValue(100,self.arg3*4)
       self._diffDataObjects(d,file)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank4(self):
       file=os.path.join(self.filebase,"tagged_function_on_boundary_rank4.nc")
       d=Data(self.arg4,FunctionOnBoundary(self.domain))
       d.setTaggedValue(1,self.arg4*2)
       d.setTaggedValue(10,self.arg4*3)
       d.setTaggedValue(100,self.arg4*4)
       self._diffDataObjects(d,file)

