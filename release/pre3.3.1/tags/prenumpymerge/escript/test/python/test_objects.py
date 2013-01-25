
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

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
          self.filename_base="."

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

import unittest
import os
import numarray
from esys.escript import *

class Test_Domain(unittest.TestCase):

   def test_getListOfTags(self): # requires self.boundary_tag_list
       tags=FunctionOnBoundary(self.domain).getListOfTags()
       self.failUnless(len(self.boundary_tag_list) == len(tags), "tag list length does not match")
       for i in self.boundary_tag_list:
           self.failUnless(i in tags, "tag %s is missing."%i)

   def test_addTags(self):
        tag1="A"
        tag2="B"
        tag3="C"
        self.domain.setTagMap(tag1,1)
        self.failUnless(self.domain.isValidTagName(tag1))
        self.failUnless(not self.domain.isValidTagName(tag2))
        self.domain.setTagMap(tag2,2)
        self.failUnless(self.domain.isValidTagName(tag1))
        self.failUnless(self.domain.isValidTagName(tag2))
        self.failUnless(not self.domain.isValidTagName(tag3))
        self.failUnless(self.domain.getTag(tag1)==1)
        self.failUnless(self.domain.getTag(tag2)==2)
        self.failUnlessRaises(RuntimeError,self.domain.getTag,tag3)

        # set tag:
        s=Scalar(0,Function(self.domain))
        r=Scalar(0,Function(self.domain))
        s.setTaggedValue(tag1,1.)
        r.setTaggedValue(1,1.)
        s.setTaggedValue(tag2,2.)
        r.setTaggedValue(2,2.)
        self.failUnlessRaises(RuntimeError,s.setTaggedValue,tag3,3.)	#tag3 does not exist
        self.failUnless(Lsup(s-r)<=0.)
        # get tag:
        names=getTagNames(self.domain)
        self.failUnless(len(names) == 6)
        self.failUnless( tag1 in names )
        self.failUnless( tag2 in names )
        self.failUnless(self.domain.isValidTagName(tag1))
        self.failUnless(self.domain.isValidTagName(tag2))
        # insert tag shortcut:
        s2=insertTaggedValues(Scalar(0,Function(self.domain)),**{ tag1 : 1., tag2 : 2.})
        self.failUnless(Lsup(s2-r)<=0.)
   def test_functionspace_ContinuousFunction(self):
        fs=ContinuousFunction(self.domain)
        self.failUnless(fs.getDomain()==self.domain)
        self.failUnless(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.failUnless(x.getFunctionSpace() == fs)
        self.failUnless(x.getShape() == (fs.getDim(),))
        self.failUnless(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.failUnless(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.failUnless(inf(x[2])>=0.)
        self.failUnless(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.failUnless(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.failUnless(sup(x[2])<=1.)

   def test_functionspace_Solution(self):
        fs=Solution(self.domain)
        self.failUnless(fs.getDomain()==self.domain)
        self.failUnless(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.failUnless(x.getFunctionSpace() == fs)
        self.failUnless(x.getShape() == (fs.getDim(),))
        self.failUnless(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.failUnless(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.failUnless(inf(x[2])>=0.)
        self.failUnless(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.failUnless(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.failUnless(sup(x[2])<=1.)

   def test_functionspace_ReducedSolution(self):
        fs=ReducedSolution(self.domain)
        self.failUnless(fs.getDomain()==self.domain)
        self.failUnless(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.failUnless(x.getFunctionSpace() == fs)
        self.failUnless(x.getShape() == (fs.getDim(),))
        self.failUnless(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.failUnless(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.failUnless(inf(x[2])>=0.)
        self.failUnless(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.failUnless(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.failUnless(sup(x[2])<=1.)

   def test_functionspace_Function(self):
        fs=Function(self.domain)
        self.failUnless(fs.getDomain()==self.domain)
        self.failUnless(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.failUnless(x.getFunctionSpace() == fs)
        self.failUnless(x.getShape() == (fs.getDim(),))
        self.failUnless(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.failUnless(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.failUnless(inf(x[2])>=0.)
        self.failUnless(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.failUnless(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.failUnless(sup(x[2])<=1.)

   def test_functionspace_ReducedFunction(self):
        fs=ReducedFunction(self.domain)
        self.failUnless(fs.getDomain()==self.domain)
        self.failUnless(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.failUnless(x.getFunctionSpace() == fs)
        self.failUnless(x.getShape() == (fs.getDim(),))
        self.failUnless(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.failUnless(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.failUnless(inf(x[2])>=0.)
        self.failUnless(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.failUnless(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.failUnless(sup(x[2])<=1.)
   def test_functionspace_FunctionOnBoundary(self):
        fs=FunctionOnBoundary(self.domain)
        self.failUnless(fs.getDomain()==self.domain)
        self.failUnless(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.failUnless(x.getFunctionSpace() == fs)
        self.failUnless(x.getShape() == (fs.getDim(),))
        self.failUnless(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.failUnless(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.failUnless(inf(x[2])>=0.)
        self.failUnless(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.failUnless(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.failUnless(sup(x[2])<=1.)

   def test_functionspace_ReducedFunctionOnBoundary(self):
        fs=ReducedFunctionOnBoundary(self.domain)
        self.failUnless(fs.getDomain()==self.domain)
        self.failUnless(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.failUnless(x.getFunctionSpace() == fs)
        self.failUnless(x.getShape() == (fs.getDim(),))
        self.failUnless(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.failUnless(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.failUnless(inf(x[2])>=0.)
        self.failUnless(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.failUnless(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.failUnless(sup(x[2])<=1.)
   #===========================================================================
class Test_SetDataPointValue(unittest.TestCase):
   arg0=9.81
   arg1=numarray.array([3.098, -3.111])
   arg2=numarray.array([[3.82, -3.81, -0.957, 0.892, -1.367], [-4.589, -1.835, -2.679, -1.517, -4.2515], [-4.909, 1.634, -2.883, -2.135, 1.187], [0.6431, 4.638, -4.616, -0.196, -4.370]])
   arg3=numarray.array([[[-2.3667, -0.040], [-4.7398, -3.2412]], [[-2.125, -2.240], [2.237, -4.279]], [[0.68720, 2.4059], [-2.4964, 3.17453]], [[-4.907, -4.9431], [-0.3604, 0.4269]], [[1.4179, 3.326], [1.356, -0.4610]], [[3.378, 2.0902], [-2.6857, 1.3585]]])
   arg4=numarray.array([[[[-3.810, -1.3597, -1.5307, 1.099], [-1.828, 0.2526, -1.4429, 2.326], [4.9732, -2.063, 1.3153, -3.809]], [[-4.8902, -4.714, 1.520, -1.931], [-3.8847, 4.3867, 1.894030, 2.432], [-1.2082, -0.8304, 2.2612, 4.6399]]], [[[-4.5922, -3.309, -0.8171, -0.7210], [2.8051, -4.93047, 0.08450, 4.3824], [0.43204, 2.1908, 4.512633, -1.8218]], [[2.2493, -4.190, -2.3893, -4.147], [-2.104, -4.635, -4.2767, -3.53151], [-2.351, -1.6614, 2.9385, 4.099]]], [[[1.710, 0.2235, -3.4917, 0.8713], [-0.2881, 4.6278, 3.603, -2.1211], [-0.565, 4.294, -2.210827, -0.37651]], [[0.6578, -2.869, -2.490, -4.789], [3.232, 2.483, 0.9531, 2.260], [-1.785, 0.42156, -1.8379, 4.212]]]])
   def test_SetDataPointValue_Function_Rank0(self):
          d=Data(self.arg0,Function(self.domain))
          d.setValueOfDataPoint(0,self.arg0*2)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg0)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg0*2)<=Lsup(self.arg0*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg0)<=Lsup(self.arg0), "wrong setting")
   def test_SetDataPointValue_Function_Rank1(self):
          d=Data(self.arg1,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,self.arg1*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_Function_Rank1_list(self):
          d=Data(self.arg1,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg1*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_Function_Rank2(self):
          d=Data(self.arg2,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,self.arg2*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_Function_Rank2_list(self):
          d=Data(self.arg2,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg2*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_Function_Rank3(self):
          d=Data(self.arg3,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,self.arg3*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_Function_Rank3_list(self):
          d=Data(self.arg3,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg3*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_Function_Rank4(self):
          d=Data(self.arg4,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,self.arg4*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")
   def test_SetDataPointValue_Function_Rank4_list(self):
          d=Data(self.arg4,Function(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg4*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")
   #===========================================================================
   def test_SetDataPointValue_ReducedFunction_Rank0(self):
          d=Data(self.arg0,ReducedFunction(self.domain))
          d.setValueOfDataPoint(0,self.arg0*2)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg0)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg0*2)<=Lsup(self.arg0*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg0)<=Lsup(self.arg0), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank1(self):
          d=Data(self.arg1,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,self.arg1*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank1_list(self):
          d=Data(self.arg1,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg1*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank2(self):
          d=Data(self.arg2,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,self.arg2*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank2_list(self):
          d=Data(self.arg2,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg2*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank3(self):
          d=Data(self.arg3,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,self.arg3*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank3_list(self):
          d=Data(self.arg3,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg3*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank4(self):
          d=Data(self.arg4,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,self.arg4*2)
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank4_list(self):
          d=Data(self.arg4,ReducedFunction(self.domain))
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.failUnlessRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,numarray.array2list(self.arg4*2))
          d_0=d.getValueOfDataPoint(0)
          d_1=d.getValueOfDataPoint(1)
          self.failUnless(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.failUnless(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")

class Test_Dump(unittest.TestCase):
   arg0=9.81
   arg1=numarray.array([3.098, -3.111])
   arg2=numarray.array([[3.82, -3.81, -0.957, 0.892, -1.367], [-4.589, -1.835, -2.679, -1.517, -4.2515], [-4.909, 1.634, -2.883, -2.135, 1.187], [0.6431, 4.638, -4.616, -0.196, -4.370]])
   arg3=numarray.array([[[-2.3667, -0.040], [-4.7398, -3.2412]], [[-2.125, -2.240], [2.237, -4.279]], [[0.68720, 2.4059], [-2.4964, 3.17453]], [[-4.907, -4.9431], [-0.3604, 0.4269]], [[1.4179, 3.326], [1.356, -0.4610]], [[3.378, 2.0902], [-2.6857, 1.3585]]])
   arg4=numarray.array([[[[-3.810, -1.3597, -1.5307, 1.099], [-1.828, 0.2526, -1.4429, 2.326], [4.9732, -2.063, 1.3153, -3.809]], [[-4.8902, -4.714, 1.520, -1.931], [-3.8847, 4.3867, 1.894030, 2.432], [-1.2082, -0.8304, 2.2612, 4.6399]]], [[[-4.5922, -3.309, -0.8171, -0.7210], [2.8051, -4.93047, 0.08450, 4.3824], [0.43204, 2.1908, 4.512633, -1.8218]], [[2.2493, -4.190, -2.3893, -4.147], [-2.104, -4.635, -4.2767, -3.53151], [-2.351, -1.6614, 2.9385, 4.099]]], [[[1.710, 0.2235, -3.4917, 0.8713], [-0.2881, 4.6278, 3.603, -2.1211], [-0.565, 4.294, -2.210827, -0.37651]], [[0.6578, -2.869, -2.490, -4.789], [3.232, 2.483, 0.9531, 2.260], [-1.785, 0.42156, -1.8379, 4.212]]]])

   def _diffDataObjects(self,d_ref,filemame, use_old_file=False):
       if not use_old_file: d_ref.dump(filemame)
       d=load(filemame, d_ref.getDomain())
       self.failUnless(not d.isEmpty(),"data in %s are empty."%filemame)
       self.failUnless(d_ref.getRank() == d.getRank(), "different rank in %s. "%filemame)
       self.failUnless(d_ref.getShape() == d.getShape(), "different shape %s. "%filemame)
       self.failUnless(d_ref.getFunctionSpace() == d.getFunctionSpace(), "wrong function space in %s."%filemame)
       self.failUnless(Lsup(d_ref-d)<=0., "different entries %s."%filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_Solution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank0.nc")
          d=Data(self.arg0,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank1.nc")
          d=Data(self.arg1,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank2.nc")
          d=Data(self.arg2,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank3.nc")
          d=Data(self.arg3,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank4.nc")
          d=Data(self.arg4,Solution(self.domain))
          self._diffDataObjects(d,filemame)
   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedSolution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank0.nc")
          d=Data(self.arg0,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank1.nc")
          d=Data(self.arg1,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank2.nc")
          d=Data(self.arg2,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank3.nc")
          d=Data(self.arg3,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank4.nc")
          d=Data(self.arg4,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
   #===========================================================================
   def test_DumpAndLoad_Constant_ContinuousFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank0.nc")
          d=Data(self.arg0,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank1.nc")
          d=Data(self.arg1,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank2.nc")
          d=Data(self.arg2,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank3.nc")
          d=Data(self.arg3,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank4.nc")
          d=Data(self.arg4,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_Function_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank0.nc")
          d=Data(self.arg0,Function(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Function_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank1.nc")
          d=Data(self.arg1,Function(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Function_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank2.nc")
          d=Data(self.arg2,Function(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Function_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank3.nc")
          d=Data(self.arg3,Function(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank0.nc")
          d=Data(self.arg0,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank1.nc")
          d=Data(self.arg1,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank2.nc")
          d=Data(self.arg2,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank3.nc")
          d=Data(self.arg3,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
   def test_DumpAndLoad_Constant_ReducedFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank4.nc")
          d=Data(self.arg4,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank0.nc")
          d=Data(self.arg0,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank1.nc")
          d=Data(self.arg1,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank2.nc")
          d=Data(self.arg2,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank3.nc")
          d=Data(self.arg3,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank4.nc")
          d=Data(self.arg4,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank0.nc")
          d=Data(self.arg0,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank1.nc")
          d=Data(self.arg1,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank2.nc")
          d=Data(self.arg2,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank3.nc")
          d=Data(self.arg3,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank4.nc")
          d=Data(self.arg4,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Solution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank0.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg0,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg0,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank1.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg1,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg1,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank2.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg2,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg2,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank3.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg3,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg3,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank4.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg4,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg4,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedSolution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank0.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg0,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg0,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank1.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg1,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg1,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank2.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg2,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg2,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank3.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg3,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg3,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank4.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg4,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg4,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank0.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg0,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg0,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank1.nc")
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg1,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg1,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank2.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg2,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg2,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank3.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg3,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg3,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank4.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg4,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg4,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Function_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank0.nc")
          d=Data(length(Function(self.domain).getX())*self.arg0,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg0,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank1.nc")
          d=Data(length(Function(self.domain).getX())*self.arg1,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg1,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank2.nc")
          d=Data(length(Function(self.domain).getX())*self.arg2,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg2,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank3.nc")
          d=Data(length(Function(self.domain).getX())*self.arg3,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg3,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank4.nc")
          d=Data(length(Function(self.domain).getX())*self.arg4,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg4,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank0.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg0,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg0,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank1.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg1,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg1,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank2.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg2,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg2,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank3.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg3,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg3,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank4.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg4,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg4,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank0.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg0,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg0,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank1.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg1,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg1,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank2.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg2,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg2,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank3.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg3,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg3,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank4.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg4,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg4,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank0.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg0,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg0,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank1.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg1,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg1,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank2.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg2,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg2,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank3.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg3,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg3,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank4.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg4,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.failUnlessRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg4,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   ## This functionspace does not currently support tags.
   ## Instead, we test that the canTag() function throws in test_canTag_Failures.
   
   #def test_DumpAndLoad_Tagged_Solution_Rank0(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank0.nc")
          #d=Data(self.arg0,Solution(self.domain))
          #d.setTaggedValue(1,self.arg0*2)
          #d.setTaggedValue(10,self.arg0*3)
          #d.setTaggedValue(100,self.arg0*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank1(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank1.nc")
          #d=Data(self.arg1,Solution(self.domain))
          #d.setTaggedValue(1,self.arg1*2)
          #d.setTaggedValue(10,self.arg1*3)
          #d.setTaggedValue(100,self.arg1*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank2(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank2.nc")
          #d=Data(self.arg2,Solution(self.domain))
          #d.setTaggedValue(1,self.arg2*2)
          #d.setTaggedValue(10,self.arg2*3)
          #d.setTaggedValue(100,self.arg2*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank3(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank3.nc")
          #d=Data(self.arg3,Solution(self.domain))
          #d.setTaggedValue(1,self.arg3*2)
          #d.setTaggedValue(10,self.arg3*3)
          #d.setTaggedValue(100,self.arg3*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank4(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank4.nc")
          #d=Data(self.arg4,Solution(self.domain))
          #d.setTaggedValue(1,self.arg4*2)
          #d.setTaggedValue(10,self.arg4*3)
          #d.setTaggedValue(100,self.arg4*4)
          #self._diffDataObjects(d,filemame)
   ##===========================================================================
   ## This functionspace does not currently support tags.
   ## Instead, we test that the canTag() function throws in test_canTag_Failures.
   
   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank0(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank0.nc")
          #d=Data(self.arg0,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg0*2)
          #d.setTaggedValue(10,self.arg0*3)
          #d.setTaggedValue(100,self.arg0*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank1(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank1.nc")
          #d=Data(self.arg1,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg1*2)
          #d.setTaggedValue(10,self.arg1*3)
          #d.setTaggedValue(100,self.arg1*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank2(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank2.nc")
          #d=Data(self.arg2,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg2*2)
          #d.setTaggedValue(10,self.arg2*3)
          #d.setTaggedValue(100,self.arg2*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank3(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank3.nc")
          #d=Data(self.arg3,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg3*2)
          #d.setTaggedValue(10,self.arg3*3)
          #d.setTaggedValue(100,self.arg3*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank4(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank4.nc")
          #d=Data(self.arg4,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg4*2)
          #d.setTaggedValue(10,self.arg4*3)
          #d.setTaggedValue(100,self.arg4*4)
          #self._diffDataObjects(d,filemame)
   ##===========================================================================
   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank0.nc")
          d=Data(self.arg0,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank1.nc")
          d=Data(self.arg1,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank2.nc")
          d=Data(self.arg2,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank3.nc")
          d=Data(self.arg3,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank4.nc")
          d=Data(self.arg4,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Tagged_Function_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank0.nc")
          d=Data(self.arg0,Function(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank1.nc")
          d=Data(self.arg1,Function(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank2.nc")
          d=Data(self.arg2,Function(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank3.nc")
          d=Data(self.arg3,Function(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank4.nc")
          d=Data(self.arg4,Function(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank0.nc")
          d=Data(self.arg0,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank1.nc")
          d=Data(self.arg1,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank2.nc")
          d=Data(self.arg2,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank3.nc")
          d=Data(self.arg3,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank4.nc")
          d=Data(self.arg4,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)
   #===========================================================================
   def test_DumpAndLoad_Tagged_ReducedFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank0.nc")
          d=Data(self.arg0,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank1.nc")
          d=Data(self.arg1,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank2.nc")
          d=Data(self.arg2,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank3.nc")
          d=Data(self.arg3,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank4.nc")
          d=Data(self.arg4,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank0.nc")
          d=Data(self.arg0,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank1.nc")
          d=Data(self.arg1,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank2.nc")
          d=Data(self.arg2,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank3.nc")
          d=Data(self.arg3,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank4.nc")
          d=Data(self.arg4,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

   def test_canTag_Failures(self):
	d=Data(self.arg0,Solution(self.domain))
	self.failUnlessRaises(RuntimeError,d.setTaggedValue,1,self.arg0*2)
	d=Data(self.arg0,ReducedSolution(self.domain))
	self.failUnlessRaises(RuntimeError,d.setTaggedValue,1,self.arg0*2)
