# $Id: run_visualization_interface.py 798 2006-08-04 01:05:36Z gross $

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

import os
import sys
import unittest
import math
import numarray
from esys.pycad import *
from esys.pycad.design import Design as Design0
from esys.pycad.gmsh import Design as GMSHDesign
from esys.pycad.Triangle import Design as TriangleDesign

try:
     PYCAD_TEST_DATA=os.environ['PYCAD_TEST_DATA']
except KeyError:
     PYCAD_TEST_DATA='.'

try:
     PYCAD_WORKDIR=os.environ['PYCAD_WORKDIR']
except KeyError:
     PYCAD_WORKDIR='.'

PYCAD_TEST_MESH_PATH=PYCAD_TEST_DATA+os.sep+"data_meshes"+os.sep
PYCAD_WORKDIR_PATH=PYCAD_WORKDIR+os.sep

def _cross(x, y):
    return numarray.array([x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]])


class Test_PyCAD_Transformations(unittest.TestCase):
   ABS_TOL=1.e-8
   def __distance(self,x,y):
       return math.sqrt(numarray.dot(x-y,x-y))
   def test_Translation_x(self):
        t=Translation([1,0,0])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([2,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([1,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([1,0,1]))<self.ABS_TOL,"s2 is wrong.")
   def test_Translation_y(self):
        t=Translation([0,1,0])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1,1,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,2,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,1,1]))<self.ABS_TOL,"s2 is wrong.")
   def test_Translation_z(self):
        t=Translation([0,0,1])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1,0,1]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,1,1]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_0_two(self):
        t=Dilation(2.)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([2,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,2,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_0_half(self):
        t=Dilation(0.5)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([0.5,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,0.5,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,0.5]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_x_two(self):
        t=Dilation(2.,[1.,0.,0.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s0_1=t([0,0,0])
        self.failUnless(isinstance(s0_1,numarray.NumArray),"s0_1 is not a numarray object.")
        self.failUnless(self.__distance(s0_1,numarray.array([-1.,0,0]))<self.ABS_TOL,"s0_1 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([-1,2,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([-1.,0,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_x_half(self):
        t=Dilation(0.5,[1.,0.,0.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s0_1=t([0,0,0])
        self.failUnless(isinstance(s0_1,numarray.NumArray),"s0_1 is not a numarray object.")
        self.failUnless(self.__distance(s0_1,numarray.array([.5,0,0]))<self.ABS_TOL,"s0_1 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.5,0.5,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.5,0,0.5]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_y_two(self):
        t=Dilation(2.,[0.,1.,0.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([2.,-1.,0]))<self.ABS_TOL,"s0 is wrong.")
        s1_1=t([0,0,0])
        self.failUnless(isinstance(s1_1,numarray.NumArray),"s1_1 is not a numarray object.")
        self.failUnless(self.__distance(s1_1,numarray.array([0.,-1.,0]))<self.ABS_TOL,"s1_1 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,1.,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,-1.,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_y_half(self):
        t=Dilation(0.5,[0.,1.,0.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([0.5,0.5,0]))<self.ABS_TOL,"s0 is wrong.")
        s1_1=t([0,0,0])
        self.failUnless(isinstance(s1_1,numarray.NumArray),"s1_1 is not a numarray object.")
        self.failUnless(self.__distance(s1_1,numarray.array([0,0.5,0]))<self.ABS_TOL,"s1_1 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,1.,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0.5,0.5]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_z_two(self):
        t=Dilation(2.,[0.,0.,1.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([2.,0.,-1.]))<self.ABS_TOL,"s0 is wrong.")
        s2_1=t([0,0,0])
        self.failUnless(isinstance(s2_1,numarray.NumArray),"s2_1 is not a numarray object.")
        self.failUnless(self.__distance(s2_1,numarray.array([0.,0.,-1.]))<self.ABS_TOL,"s2_1 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,2.,-1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0.,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_z_half(self):
        t=Dilation(0.5,[0.,0.,1.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([0.5,0.,0.5]))<self.ABS_TOL,"s0 is wrong.")
        s2_1=t([0,0,0])
        self.failUnless(isinstance(s2_1,numarray.NumArray),"s2_1 is not a numarray object.")
        self.failUnless(self.__distance(s2_1,numarray.array([0,0,0.5]))<self.ABS_TOL,"s2_1 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,0.5,0.5]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0.,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Reflection_x_offset0(self):
        t=Reflection([1.,0.,0.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([-1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([-1.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_x_offset2(self):
        t=Reflection([-2.,0.,0.],offset=-4)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([3.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([4,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([4,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([3.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_x_offset2_vector(self):
        t=Reflection([1.,0.,0.],offset=[2,0,0])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([3.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([4,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([4,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([3.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_y_offset0(self):
        t=Reflection([0.,1.,0.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,-1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([1.,-2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_y_offset2(self):
        t=Reflection([0.,-2.,0.],offset=-4)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,4,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,3,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,4,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([1.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_y_offset2_vector(self):
        t=Reflection([0.,1.,0.],offset=[0,2,0])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,4,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,3,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,4,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([1.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_z_offset0(self):
        t=Reflection([0.,0.,1.])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,-1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([1.,2,-3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_z_offset2(self):
        t=Reflection([0.,0.,-2.],offset=-4)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,4.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,1,4]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,3]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([1.,2,1]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_z_offset2_vector(self):
        t=Reflection([0.,0.,1.],offset=[0,0,2])
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,4.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0,1,4]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0,0,3]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.failUnless(isinstance(s,numarray.NumArray),"s is not a numarray object.")
        self.failUnless(self.__distance(s,numarray.array([1.,2,1]))<self.ABS_TOL,"s is wrong.")
   def test_Rotatation_x_90_0(self):
        t=Rotatation(axis=[1.,0.,0.],point=[1.,0.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,0,1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,-1.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_x_30_0(self):
        t=Rotatation(axis=[1.,0.,0.],point=[1.,0.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,0]),numarray.array([1.,0.,0.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([1.,0.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_x_330_0(self):
        t=Rotatation(axis=[1.,0.,0.],point=[1.,0.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]-math.cos(330*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,0]),numarray.array([1.,0.,0.]))>0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([1.,0.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_x_90(self):
        t=Rotatation(axis=[-1.,0.,0.],point=[2.,0.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,0,-1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,1.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_x_30(self):
        t=Rotatation(axis=[-1.,0.,0.],point=[1.,0.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,0]),numarray.array([-1.,0.,0.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([-1.,0.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_x_330(self):
        t=Rotatation(axis=[-1.,0.,0.],point=[1.,0.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]-math.cos(330*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,0]),numarray.array([-1.,0.,0.]))>0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([-1.,0.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_y_90_0(self):
        t=Rotatation(axis=[0.,1.,0.],point=[0.,1.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([0.,0,-1.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,5,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([1,0.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_y_30_0(self):
        t=Rotatation(axis=[0.,1.,0.],point=[0.,1.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,1.,0.]))<0.,"s0 has wrong orientation.")
        s1=t([0,5,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,5,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([0.,1.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_y_330_0(self):
        t=Rotatation(axis=[0.,1.,0.],point=[0.,1.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,1.,0.]))>0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([0.,1.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_y_90(self):
        t=Rotatation(axis=[0.,-1.,0.],point=[0.,2.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([0.,0,1.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,5,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([-1,0.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_y_30(self):
        t=Rotatation(axis=[0.,-1.,0.],point=[0.,2.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,-1.,0.]))<0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([0.,-1.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_y_330(self):
        t=Rotatation(axis=[0.,-1.,0.],point=[0.,2.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,-1.,0.]))>0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[0,0,1]),numarray.array([0.,-1.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_z_90_0(self):
        t=Rotatation(axis=[0.,0.,1.],point=[0.,0.,1.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([0.,1,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([-5.,0,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_30_0(self):
        t=Rotatation(axis=[0.,0.,1.],point=[0.,0.,1.],angle=30*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,0.,1.]))<0.,"s0 has wrong orientation.")
        s1=t([0,5,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-5.**2)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]/5.-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,5,0]),numarray.array([0.,0.,1.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_330_0(self):
        t=Rotatation(axis=[0.,0.,1.],point=[0.,0.,1.],angle=330*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,0.,1.]))>0.,"s0 has wrong orientation.")
        s1=t([0,5,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-5.**2)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]/5.-math.cos(330*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,0]),numarray.array([0.,0.,1.]))>0.,"s1 has wrong orientation.")
   def test_Rotatation_z_90(self):
        t=Rotatation(axis=[0.,0.,-1.],point=[0.,0.,2.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([0.,-1,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([5.,0,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_30(self):
        t=Rotatation(axis=[0.,0.,-1.],point=[0.,0.,2.],angle=30*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,0.,-1.]))<0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,0]),numarray.array([0.,0.,-1.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_330(self):
        t=Rotatation(axis=[0.,0.,-1.],point=[0.,0.,2.],angle=330*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,0,0]),numarray.array([0.,0.,-1.]))>0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,0]),numarray.array([0.,0.,-1.]))>0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_x_90_1(self):
        t=Rotatation(point=[0.,0.,1.],axis=[1.,0.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,1,1.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,1,2.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_y_90_1(self):
        t=Rotatation(point=[1.,0.,0.],axis=[0.,1.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([1.,1,1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([2.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_90_1(self):
        t=Rotatation(point=[0.,1.,0.],axis=[0.,0.,1.],angle=90*DEG)
        s0=t([1,0,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(self.__distance(s0,numarray.array([1.,2,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(self.__distance(s1,numarray.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(self.__distance(s2,numarray.array([1.,1,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_diag_90_0(self):
        t=Rotatation(axis=[1.,1.,1.],angle=90*DEG)
        s0=t([1,-1,0])
        self.failUnless(isinstance(s0,numarray.NumArray),"s0 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s0,s0)-2.)<self.ABS_TOL,"s0 length is wrong.")
        self.failUnless(abs(numarray.dot(s0,numarray.array([1,-1,0])))<self.ABS_TOL,"s0 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s0,[1,-1,0]),numarray.array([1.,1.,1.]))<0.,"s0 has wrong orientation.")
        s1=t([0,1,-1])
        self.failUnless(isinstance(s1,numarray.NumArray),"s1 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s1,s1)-2.)<self.ABS_TOL,"s1 length is wrong.")
        self.failUnless(abs(numarray.dot(s1,numarray.array([0,1,-1])))<self.ABS_TOL,"s1 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s1,[0,1,-1]),numarray.array([1.,1.,1.]))<0.,"s1 has wrong orientation.")
        s2=t([-1,0,1])
        self.failUnless(isinstance(s2,numarray.NumArray),"s2 is not a numarray object.")
        self.failUnless(abs(numarray.dot(s2,s2)-2.)<self.ABS_TOL,"s2 length is wrong.")
        self.failUnless(abs(numarray.dot(s2,numarray.array([-1,0,1])))<self.ABS_TOL,"s2 angle is wrong.")
        self.failUnless(numarray.dot(_cross(s2,[-1,0,1]),numarray.array([1.,1.,1.]))<0.,"s2 has wrong orientation.")
        s3=t([1,1,1])
        self.failUnless(isinstance(s3,numarray.NumArray),"s3 is not a numarray object.")
        self.failUnless(self.__distance(s3,numarray.array([1.,1,1.]))<self.ABS_TOL,"s3 is wrong.")

class Test_PyCAD_Primitives(unittest.TestCase):
   def setUp(self):
         resetGlobalPrimitiveIdCounter()

   def test_Primitive(self):
         p=Primitive()

         id=p.getID()
         self.failUnless(isinstance(id,int),"id number is not an integer")
         self.failUnless(not id==Primitive().getID(),"id number is not unique")

         self.failUnless(p==p.getUnderlyingPrimitive(),"getUnderlyingPrimitive does not return self.")

   def test_ReversePrimitive(self):
         p=Primitive()
      
         rp=ReversePrimitive(p)
         self.failUnless(p.getID()==rp.getID(),"reverse primitive does not have same id like source")
         self.failUnless(p==rp.getUnderlyingPrimitive(),"getUnderlyingPrimitive does return source.")
         self.failUnless(p == -rp,"reverse or reverse does not return source.")
           
   def test_Point(self):
       p=Point(1.,2.,3.,local_scale=9.)
       
       id=p.getID()
       self.failUnless(isinstance(id,int),"id number is not an integer")
       self.failUnless(not id==Primitive().getID(),"id number is not unique")
           
       # check reverse point 
       self.failUnless(p == -p,"reverse is not working.")
       
       # check history:
       hs=p.getPrimitives()
       self.failUnless(len(hs)==1,"history must have length 1.")
       self.failUnless(p in hs,"history must contain point p")

       # check incolved points:
       ps=p.getConstructionPoints()
       self.failUnless(len(ps)==1,"point set must have length 1.")
       self.failUnless(p in ps,"point set must contain point p")

       # check coordinates:
       c=p.getCoordinates()
       self.failUnless(isinstance(c,numarray.NumArray),"coordinates are not a numarray object.")
       self.failUnless(c[0]==1.,"x coordinate is not 1.")
       self.failUnless(c[1]==2.,"y coordinate is not 2.")
       self.failUnless(c[2]==3.,"z coordinate is not 3.")
 
       # reset coordinates:
       p.setCoordinates([-1.,-2.,-3.])
       c=p.getCoordinates()
       self.failUnless(isinstance(c,numarray.NumArray),"new coordinates are not a numarray object.")
       self.failUnless(c[0]==-1.,"new x coordinate is not -1.")
       self.failUnless(c[1]==-2.,"new y coordinate is not -2.")
       self.failUnless(c[2]==-3.,"new z coordinate is not -3.")

       # check for a colocated point:
       self.failUnless(p.isColocated(Point(-1.,-2.,-3.)),"colocation not detected.")
       self.failUnless(not p.isColocated(numarray.array([-1.,-2.,-3.])),"colocation with numarray representation not detected.")
       self.failUnless(not p.isColocated(Point(1.,-2.,-3.)),"false colocation detected.")
       self.failUnless(not p.isColocated(Point(0.,0.,0.)),"false colocation with origin detected.")

       # check for local length scale
       l=p.getLocalScale()
       self.failUnless(l==9.,"refinement scale is not 9.")

       # check for new local length scale
       p.setLocalScale(3.)
       l=p.getLocalScale()
       self.failUnless(l==3.,"new refinement scale is not 3.")

       # negative value shouldn't work.
       self.failUnlessRaises(ValueError,p.setLocalScale,-3.)

       # copy:
       an_other_p=p.copy()
       self.failUnless(isinstance(an_other_p ,Point),"copy is not a point")
       self.failUnless(not an_other_p.getID() == p.getID(),"copy has same Id")
       self.failUnless(p.isColocated(an_other_p),"p is not colocated with its copy.")
       self.failUnless(an_other_p.isColocated(p),"the copy is not colocated with p.")
       self.failUnless(an_other_p.getLocalScale()==3.,"copy has wrong local scale.")
      
       # modify by Transformation:
       p.modifyBy(Dilation(-1))
       self.failUnless(p.isColocated(Point(1.,2.,3.)),"in-place transformation failed")
       
       # apply Transformation:
       dil_p=p.apply(Dilation(4))
       self.failUnless(dil_p.isColocated(Point(4.,8.,12.)),"applying transformation failed")
       self.failUnless(not dil_p.getID() == p.getID(),"transformed point has same Id")
       self.failUnless(dil_p.getLocalScale()==3.,"transformed point  has wrong local scale.")
        
       # overloaded add:
       shift_p=p+[1,1,1]
       self.failUnless(shift_p.isColocated(Point(2,3.,4)),"applying shift by list failed")
       self.failUnless(not shift_p.getID() == p.getID(),"shift by list has same Id")
       self.failUnless(shift_p.getLocalScale()==3.,"shift by list has wrong local scale.")

       shift_p=p+numarray.array([1,1,1])
       self.failUnless(shift_p.isColocated(Point(2,3.,4)),"applying shift by numarray failed")
       self.failUnless(not shift_p.getID() == p.getID(),"shift by numarray has same Id")
       self.failUnless(shift_p.getLocalScale()==3.,"shift by numarray has wrong local scale.")
       # overloaded minus
       shift_p=p-[1,1,1]
       self.failUnless(shift_p.isColocated(Point(0,1,2.)),"applying shift by -list failed")
       self.failUnless(not shift_p.getID() == p.getID(),"shift by -list has same Id")
       self.failUnless(shift_p.getLocalScale()==3.,"shift by -list has wrong local scale.")

       shift_p=p-numarray.array([1,1,1])
       self.failUnless(shift_p.isColocated(Point(0,1,2.)),"applying shift by -numarray failed")
       self.failUnless(not shift_p.getID() == p.getID(),"shift by -numarray has same Id")
       self.failUnless(shift_p.getLocalScale()==3.,"shift by -numarray has wrong local scale.")
       # overloaded inplace add:
       p+=[1,1,1]
       self.failUnless(p.isColocated(Point(2,3.,4)),"modification by list shift failed")

       p+=numarray.array([1,1,1])
       self.failUnless(p.isColocated(Point(3,4,5)),"modification by numarray shift failed")

       # overloaded inplace add:
       p-=[1,1,1]
       self.failUnless(p.isColocated(Point(2,3,4)),"modification by -list shift failed")

       p-=numarray.array([1,1,1])
       self.failUnless(p.isColocated(Point(1,2.,3)),"modification by -numarray shift failed")

       #overloaded multiplication:
       mult_p=2*p
       self.failUnless(mult_p.isColocated(Point(2,4,6)),"applying int factor failed")
       self.failUnless(not mult_p.getID() == p.getID(),"shift by int factor has same Id")
       self.failUnless(mult_p.getLocalScale()==3.,"shift by int factor has wrong local scale.")

       mult_p=2.*p
       self.failUnless(mult_p.isColocated(Point(2,4,6)),"applying float factor failed")
       self.failUnless(not mult_p.getID() == p.getID(),"shift by float factor has same Id")
       self.failUnless(mult_p.getLocalScale()==3.,"shift by float factor has wrong local scale.")

       mult_p=Dilation(2)*p
       self.failUnless(mult_p.isColocated(Point(2,4,6)),"applying Dilation factor failed")
       self.failUnless(not mult_p.getID() == p.getID(),"shift by Dilation factor has same Id")
       self.failUnless(mult_p.getLocalScale()==3.,"shift by Dilation factor has wrong local scale.")

       #overloaded inplace multiplication:
       p*=2
       self.failUnless(p.isColocated(Point(2,4,6)),"applying in-place int factor failed")

       p*=2.
       self.failUnless(p.isColocated(Point(4,8,12)),"applying in-place float factor failed")

       p*=Dilation(2)
       self.failUnless(p.isColocated(Point(8,16,24)),"applying in-place Dilation factor failed")

   def test_Spline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
 
        self.failUnlessRaises(ValueError,Spline,p0)
        c=Spline(p0,p1,p2,p3)

        self.failUnless(len(c) == 4, "wrong spline curve length")
        self.failUnless(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.failUnless(c.getEndPoint()==p3, "wrong end point of spline curve")

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p1),"spline is colocated with point.")
        self.failUnless(not c.isColocated(Spline(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.failUnless(not c.isColocated(Spline(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.failUnless(c.isColocated(Spline(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(Spline(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.failUnless(co[0]==p0, "1st control point is wrong.")
        self.failUnless(co[1]==p1, "2nd control point is wrong.")
        self.failUnless(co[2]==p2, "3rd control point is wrong.")
        self.failUnless(co[3]==p3, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.failUnless(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.failUnless(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.failUnless(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.failUnless(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 5, "number of primitives in history is wrong.")
        self.failUnless(p0 in h, "missing p0 in history.")
        self.failUnless(p1 in h, "missing p1 in history.")
        self.failUnless(p2 in h, "missing p2 in history.")
        self.failUnless(p3 in h, "missing p3 in history.")
        self.failUnless(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.failUnless(not cp == c, "copy returns same spline curve.")
        self.failUnless(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.failUnless(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.failUnless(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.failUnless(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.failUnless(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.failUnless(c.isColocated(Spline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.failUnless(p0 == cp[0],"1st new point after Dilation.")
        self.failUnless(p1 == cp[1],"2nd new point after Dilation.")
        self.failUnless(p2 == cp[2],"3rd new point after Dilation.")
        self.failUnless(p3 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.failUnless(dc.isColocated(Spline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.failUnless(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.failUnless(dccp[0].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.failUnless(dccp[1].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.failUnless(dccp[2].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p3 == dccp[3],"4th point of Dilation is identical to source.")
        self.failUnless(dccp[3].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")
   
   def test_ReverseSpline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
 
        CC0=Spline(p0,p1,p2,p3)
        c=-CC0

        self.failUnless(len(c) == 4, "wrong reverse spline curve length")
        self.failUnless(c.getStartPoint()==p3, "wrong start point of reverse spline curve")
        self.failUnless(c.getEndPoint()==p0, "wrong end point of reverse spline curve")

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p1),"reverse spline is colocated with point.")
        self.failUnless(not c.isColocated(Spline(p0,p1,p2)),"reverse spline is colocated with spline of different length.")
        self.failUnless(not c.isColocated(Spline(p0,p1,p4,p3)),"reverse spline is colocated with spline with different point.")
        self.failUnless(c.isColocated(Spline(p0,p1,p2,p3)),"reverse spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(c.isColocated(Spline(p3,p2,p1,p0)),"reverse spline is not colocated with spline with same points.")
        self.failUnless(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.failUnless(co[0]==p3, "1st control point is wrong.")
        self.failUnless(co[1]==p2, "2nd control point is wrong.")
        self.failUnless(co[2]==p1, "3rd control point is wrong.")
        self.failUnless(co[3]==p0, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.failUnless(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.failUnless(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.failUnless(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.failUnless(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 5, "number of primitives in history is wrong.")
        self.failUnless(p0 in h, "missing p0 in history.")
        self.failUnless(p1 in h, "missing p1 in history.")
        self.failUnless(p2 in h, "missing p2 in history.")
        self.failUnless(p3 in h, "missing p3 in history.")
        self.failUnless(CC0 in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.failUnless(not cp == c, "copy returns same spline curve.")
        self.failUnless(not cp == CC0, "copy returns same spline curve.")
        self.failUnless(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.failUnless(not p3 == cpcp[0],"1st point of deep copy and souce are the same.")
        self.failUnless(not p2 == cpcp[1],"2st point of deep copy and source are the same.")
        self.failUnless(not p1 == cpcp[2],"3st point of deep copy and source are the same.")
        self.failUnless(not p0 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.failUnless(c.isColocated(Spline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.failUnless(p3 == cp[0],"1st new point after Dilation.")
        self.failUnless(p2 == cp[1],"2nd new point after Dilation.")
        self.failUnless(p1 == cp[2],"3rd new point after Dilation.")
        self.failUnless(p0 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.failUnless(dc.isColocated(Spline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.failUnless(dccp[0].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")
        self.failUnless(dccp[1].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.failUnless(dccp[2].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.failUnless(dccp[3].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")

   def test_BezierCurve(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
 
        self.failUnlessRaises(ValueError,BezierCurve,p0)
        c=BezierCurve(p0,p1,p2,p3)

        self.failUnless(len(c) == 4, "wrong spline curve length")
        self.failUnless(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.failUnless(c.getEndPoint()==p3, "wrong end point of spline curve")

        self.failUnless(not c.isColocated(p1),"spline is colocated with point.")
        self.failUnless(not c.isColocated(BezierCurve(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.failUnless(not c.isColocated(BezierCurve(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.failUnless(c.isColocated(BezierCurve(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(BezierCurve(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.failUnless(co[0]==p0, "1st control point is wrong.")
        self.failUnless(co[1]==p1, "2nd control point is wrong.")
        self.failUnless(co[2]==p2, "3rd control point is wrong.")
        self.failUnless(co[3]==p3, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.failUnless(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.failUnless(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.failUnless(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.failUnless(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 5, "number of primitives in history is wrong.")
        self.failUnless(p0 in h, "missing p0 in history.")
        self.failUnless(p1 in h, "missing p1 in history.")
        self.failUnless(p2 in h, "missing p2 in history.")
        self.failUnless(p3 in h, "missing p3 in history.")
        self.failUnless(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.failUnless(not cp == c, "copy returns same spline curve.")
        self.failUnless(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.failUnless(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.failUnless(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.failUnless(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.failUnless(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.failUnless(c.isColocated(BezierCurve(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.failUnless(p0 == cp[0],"1st new point after Dilation.")
        self.failUnless(p1 == cp[1],"2nd new point after Dilation.")
        self.failUnless(p2 == cp[2],"3rd new point after Dilation.")
        self.failUnless(p3 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.failUnless(dc.isColocated(BezierCurve(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.failUnless(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.failUnless(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.failUnless(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.failUnless(not p3 == dccp[3],"4th point of Dilation is identical to source.")

   def test_BSpline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
 
        self.failUnlessRaises(ValueError,BSpline,p0)
        c=BSpline(p0,p1,p2,p3)

        self.failUnless(len(c) == 4, "wrong spline curve length")
        self.failUnless(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.failUnless(c.getEndPoint()==p3, "wrong end point of spline curve")

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p1),"spline is colocated with point.")
        self.failUnless(not c.isColocated(BSpline(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.failUnless(not c.isColocated(BSpline(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.failUnless(c.isColocated(BSpline(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(BSpline(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.failUnless(co[0]==p0, "1st control point is wrong.")
        self.failUnless(co[1]==p1, "2nd control point is wrong.")
        self.failUnless(co[2]==p2, "3rd control point is wrong.")
        self.failUnless(co[3]==p3, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.failUnless(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.failUnless(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.failUnless(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.failUnless(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 5, "number of primitives in history is wrong.")
        self.failUnless(p0 in h, "missing p0 in history.")
        self.failUnless(p1 in h, "missing p1 in history.")
        self.failUnless(p2 in h, "missing p2 in history.")
        self.failUnless(p3 in h, "missing p3 in history.")
        self.failUnless(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.failUnless(not cp == c, "copy returns same spline curve.")
        self.failUnless(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.failUnless(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.failUnless(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.failUnless(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.failUnless(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.failUnless(c.isColocated(BSpline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.failUnless(p0 == cp[0],"1st new point after Dilation.")
        self.failUnless(p1 == cp[1],"2nd new point after Dilation.")
        self.failUnless(p2 == cp[2],"3rd new point after Dilation.")
        self.failUnless(p3 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.failUnless(dc.isColocated(BSpline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.failUnless(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.failUnless(dccp[0].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.failUnless(dccp[1].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.failUnless(dccp[2].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p3 == dccp[3],"4th point of Dilation is identical to source.")
        self.failUnless(dccp[3].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")

   def test_ReverseBSpline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
 
        CC0=BSpline(p0,p1,p2,p3)
        c=-CC0

        self.failUnless(len(c) == 4, "wrong spline curve length")
        self.failUnless(c.getStartPoint()==p3, "wrong start point of spline curve")
        self.failUnless(c.getEndPoint()==p0, "wrong end point of spline curve")

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p1),"spline is colocated with point.")
        self.failUnless(not c.isColocated(BSpline(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.failUnless(not c.isColocated(BSpline(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.failUnless(c.isColocated(BSpline(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(BSpline(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.failUnless(co[0]==p3, "1st control point is wrong.")
        self.failUnless(co[1]==p2, "2nd control point is wrong.")
        self.failUnless(co[2]==p1, "3rd control point is wrong.")
        self.failUnless(co[3]==p0, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.failUnless(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.failUnless(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.failUnless(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.failUnless(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 5, "number of primitives in history is wrong.")
        self.failUnless(p0 in h, "missing p0 in history.")
        self.failUnless(p1 in h, "missing p1 in history.")
        self.failUnless(p2 in h, "missing p2 in history.")
        self.failUnless(p3 in h, "missing p3 in history.")
        self.failUnless(CC0 in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.failUnless(not cp == c, "copy returns same spline curve.")
        self.failUnless(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.failUnless(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.failUnless(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.failUnless(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.failUnless(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.failUnless(c.isColocated(BSpline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.failUnless(p3 == cp[0],"1st new point after Dilation.")
        self.failUnless(p2 == cp[1],"2nd new point after Dilation.")
        self.failUnless(p1 == cp[2],"3rd new point after Dilation.")
        self.failUnless(p0 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.failUnless(dc.isColocated(BSpline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.failUnless(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.failUnless(dccp[0].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.failUnless(dccp[1].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.failUnless(dccp[2].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p3 == dccp[3],"4th point of Dilation is identical to source.")
        self.failUnless(dccp[3].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")

   def test_LineSegment(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p4=Point(1,2,3)
 
        self.failUnlessRaises(TypeError,Line,p0)
        self.failUnlessRaises(TypeError,Line,p0,p1,p4)

        c=Line(p0,p1)

        self.failUnless(len(c) == 2, "wrong spline curve length")
        self.failUnless(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.failUnless(c.getEndPoint()==p1, "wrong end point of spline curve")

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p1),"spline is colocated with point.")
        self.failUnless(not c.isColocated(Line(p0,p4)),"spline is colocated with spline with different point.")
        self.failUnless(c.isColocated(Line(p0,p1)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(Line(p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(p0,p1,p4)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.failUnless(co[0]==p0, "1st control point is wrong.")
        self.failUnless(co[1]==p1, "2nd control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.failUnless(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.failUnless(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 3, "number of primitives in history is wrong.")
        self.failUnless(p0 in h, "missing p0 in history.")
        self.failUnless(p1 in h, "missing p1 in history.")
        self.failUnless(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.failUnless(not cp == c, "copy returns same spline curve.")
        self.failUnless(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.failUnless(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.failUnless(not p1 == cpcp[1],"2st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.failUnless(c.isColocated(Line(Point(0,0,0),Point(-1,-1,-1))),"inplace dilation is wrong.")
        self.failUnless(p0 == cp[0],"1st new point after Dilation.")
        self.failUnless(p1 == cp[1],"2nd new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.failUnless(dc.isColocated(Line(Point(0,0,0),Point(1,1,1))),"dilation is wrong.")
        self.failUnless(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.failUnless(dccp[0].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.failUnless(dccp[1].isColocated(Point(1,1,1)),"2st point of Dilation is is wrongly located.")

   def test_ReverseLineSegment(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p4=Point(1,2,3)
 
        self.failUnlessRaises(TypeError,Line,p0)
        self.failUnlessRaises(TypeError,Line,p0,p1,p4)

        CC0=Line(p0,p1)
        c=-CC0

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(len(c) == 2, "wrong spline curve length")
        self.failUnless(c.getStartPoint()==p1, "wrong start point of spline curve")
        self.failUnless(c.getEndPoint()==p0, "wrong end point of spline curve")

        self.failUnless(not c.isColocated(p1),"spline is colocated with point.")
        self.failUnless(not c.isColocated(Line(p0,p4)),"spline is colocated with spline with different point.")
        self.failUnless(c.isColocated(Line(p0,p1)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(Line(p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(p0,p1,p4)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.failUnless(co[0]==p1, "1st control point is wrong.")
        self.failUnless(co[1]==p0, "2nd control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.failUnless(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.failUnless(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 3, "number of primitives in history is wrong.")
        self.failUnless(p0 in h, "missing p0 in history.")
        self.failUnless(p1 in h, "missing p1 in history.")
        self.failUnless(CC0 in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.failUnless(not cp == c, "copy returns same spline curve.")
        self.failUnless(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.failUnless(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.failUnless(not p1 == cpcp[1],"2st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.failUnless(c.isColocated(Line(Point(0,0,0),Point(-1,-1,-1))),"inplace dilation is wrong.")
        self.failUnless(p1 == cp[0],"1st new point after Dilation.")
        self.failUnless(p0 == cp[1],"2nd new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.failUnless(dc.isColocated(Line(Point(0,0,0),Point(1,1,1))),"dilation is wrong.")
        self.failUnless(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.failUnless(dccp[0].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.failUnless(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.failUnless(dccp[1].isColocated(Point(0,0,0)),"2st point of Dilation is is wrongly located.")

   def test_Arc(self):
        center=Point(0,0,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
        p4=Point(10,2,3)
 
        self.failUnlessRaises(TypeError,Arc,Primitive())

        c=Arc(center,p_start,p_end)

        self.failUnless(c.getCenterPoint()==center, "wrong center point")
        self.failUnless(c.getStartPoint()==p_start, "wrong start point")
        self.failUnless(c.getEndPoint()==p_end, "wrong end point")

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p4),"spline is colocated with point.")
        self.failUnless(not c.isColocated(Arc(p4,p_start,p_end)),"spline is colocated with spline with differnt center point.")
        self.failUnless(not c.isColocated(Arc(center,p4,p_end)),"spline is colocated with spline with differnt start point.")
        self.failUnless(not c.isColocated(Arc(center,p_start,p4)),"spline is colocated with spline with differnt end point.")
        self.failUnless(c.isColocated(Arc(center,p_start,p_end)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(Arc(center,p_end,p_start)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(center,p_start,p_end)),"spline curve is identified with curve.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 4, "number of primitives in history is wrong.")
        self.failUnless(center in h, "missing center in history.")
        self.failUnless(p_start in h, "missing p_start in history.")
        self.failUnless(p_end in h, "missing p_end in history.")
        self.failUnless(c in h, "missing spline curve in history.")


        c.setLocalScale(3.)
        self.failUnless(c.getCenterPoint().getLocalScale() == 3., "new local scale of center point is wrong.")
        self.failUnless(c.getStartPoint().getLocalScale() == 3., "new local scale of start point is wrong.")
        self.failUnless(c.getEndPoint().getLocalScale() == 3., "new local scale of end point is wrong.")

        cp=c.copy()
        self.failUnless(isinstance(cp,Arc), "copy returns is not an arc.")
        self.failUnless(not cp == c, "copy returns same arc.")
        self.failUnless(cp.isColocated(Arc(center,p_start,p_end)),"arc is not colocated with its copy.")
        self.failUnless(not cp.getCenterPoint()==center, "deep copy has same center point like source")
        self.failUnless(not cp.getStartPoint()==p_start, "deep copy has same start point like source")
        self.failUnless(not cp.getEndPoint()==p_end, "deep copy has same end point like source")

        c.modifyBy(Dilation(-1.))
        self.failUnless(c.isColocated(Arc(Point(0,0,0),Point(-1,-1,-1),Point(-1,-2,-3))),"inplace dilation is wrong.")
        self.failUnless(c.getCenterPoint() == center,"wrong center point after dilation.")
        self.failUnless(c.getStartPoint() == p_start,"wrong start point after dilation.")
        self.failUnless(c.getEndPoint() == p_end,"wrong end point after dilation.")

        dc=c.apply(Dilation(-1.))
        self.failUnless(dc.isColocated(Arc(Point(0,0,0),Point(1,1,1),Point(1,2,3))),"dilation is wrong.")
        self.failUnless(not dc.getCenterPoint() == center,"center point of dilation is identical to source.")
        self.failUnless(dc.getCenterPoint().isColocated(Point(0,0,0)),"center point of dilation is wrong.")
        self.failUnless(not dc.getStartPoint() == p_start,"start point of dilation is identical to source.")
        self.failUnless(dc.getStartPoint().isColocated(Point(1,1,1)),"start point of dilation is wrong.")
        self.failUnless(not dc.getEndPoint() == p_end,"end point of dilation is identical to source.")
        self.failUnless(dc.getEndPoint().isColocated(Point(1,2,3)),"end point of dilation is wrong.")

   def test_ReverseArc(self):
        center=Point(0,0,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
        p4=Point(10,2,3)
 
        self.failUnlessRaises(TypeError,Arc,Primitive())

        CC0=Arc(center,p_start,p_end)
        c=-CC0

        self.failUnless(c.getCenterPoint()==center, "wrong center point")
        self.failUnless(c.getStartPoint()==p_end, "wrong start point")
        self.failUnless(c.getEndPoint()==p_start, "wrong end point")

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p4),"spline is colocated with point.")
        self.failUnless(not c.isColocated(Arc(p4,p_start,p_end)),"spline is colocated with spline with differnt center point.")
        self.failUnless(not c.isColocated(Arc(center,p4,p_end)),"spline is colocated with spline with differnt start point.")
        self.failUnless(not c.isColocated(Arc(center,p_start,p4)),"spline is colocated with spline with differnt end point.")
        self.failUnless(c.isColocated(Arc(center,p_start,p_end)),"spline is not colocated with spline with same points.")
        self.failUnless(c.isColocated(Arc(center,p_end,p_start)),"spline is not colocated with spline with same points but opposite direction.")
        self.failUnless(not c.isColocated(Curve(center,p_start,p_end)),"spline curve is identified with curve.")

        h=c.getPrimitives()
        self.failUnless(len(h) == 4, "number of primitives in history is wrong.")
        self.failUnless(center in h, "missing center in history.")
        self.failUnless(p_start in h, "missing p_start in history.")
        self.failUnless(p_end in h, "missing p_end in history.")
        self.failUnless(CC0 in h, "missing spline curve in history.")


        c.setLocalScale(3.)
        self.failUnless(c.getCenterPoint().getLocalScale() == 3., "new local scale of center point is wrong.")
        self.failUnless(c.getStartPoint().getLocalScale() == 3., "new local scale of start point is wrong.")
        self.failUnless(c.getEndPoint().getLocalScale() == 3., "new local scale of end point is wrong.")

        cp=c.copy()
        self.failUnless(isinstance(cp,ReverseArc), "copy returns is not an arc.")
        self.failUnless(not cp == c, "copy returns same arc.")
        self.failUnless(cp.isColocated(Arc(center,p_end,p_start)),"arc is not colocated with its copy.")
        self.failUnless(not cp.getCenterPoint()==center, "deep copy has same center point like source")
        self.failUnless(not cp.getStartPoint()==p_start, "deep copy has same start point like source")
        self.failUnless(not cp.getEndPoint()==p_end, "deep copy has same end point like source")

        c.modifyBy(Dilation(-1.))
        self.failUnless(c.isColocated(Arc(Point(0,0,0),Point(-1,-1,-1),Point(-1,-2,-3))),"inplace dilation is wrong.")
        self.failUnless(c.getCenterPoint() == center,"wrong center point after dilation.")
        self.failUnless(c.getStartPoint() == p_end,"wrong start point after dilation.")
        self.failUnless(c.getEndPoint() == p_start,"wrong end point after dilation.")

        dc=c.apply(Dilation(-1.))
        self.failUnless(dc.isColocated(Arc(Point(0,0,0),Point(1,1,1),Point(1,2,3))),"dilation is wrong.")
        self.failUnless(not dc.getCenterPoint() == center,"center point of dilation is identical to source.")
        self.failUnless(dc.getCenterPoint().isColocated(Point(0,0,0)),"center point of dilation is wrong.")
        self.failUnless(not dc.getStartPoint() == p_start,"start point of dilation is identical to source.")
        self.failUnless(dc.getStartPoint().isColocated(Point(1,2,3)),"start point of dilation is wrong.")
        self.failUnless(not dc.getEndPoint() == p_end,"end point of dilation is identical to source.")
        self.failUnless(dc.getEndPoint().isColocated(Point(1,1,1)),"end point of dilation is wrong.")

   def test_CurveLoop(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        p5=Point(10,20,3)
        p6=Point(1,2,30)

        l01=Line(p0,p1)
        l12=Arc(p3,p1,p2)
        l20=Spline(p2,p4,p0)

        lx=Line(p2,p3)
        ly=Line(p3,p1)

        c=CurveLoop(l01,l12,l20)
        # self.failUnlessRaises(ValueError,CurveLoop,l01,lx,l20)
        # self.failUnlessRaises(ValueError,CurveLoop,l01,l20,l20)
        # self.failUnlessRaises(ValueError,CurveLoop,l01,l20,ly)

        c=CurveLoop(l01,l20,l12)
        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p4),"CurveLoop is colocated with point.")
        self.failUnless(c.isColocated(c),"CurveLoop is not colocated with its self.")
        self.failUnless(c.isColocated(CurveLoop(l01,l12,l20)),"CurveLoop is not colocated with its copy.")
        self.failUnless(c.isColocated(CurveLoop(l20,l01,l12)),"CurveLoop is not colocated with its copy with shifted points.")
        self.failUnless(c.isColocated(CurveLoop(l20,l12,l01)),"CurveLoop is not colocated with its copy with shuffled points.")
        self.failUnless(not c.isColocated(CurveLoop(lx,ly,l12)),"CurveLoop is colocated with different CurveLoop.")

        self.failUnless(len(c) == 3, "wrong length")

        c.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")


        cc=c.getCurves()
        self.failUnless(len(cc) == 3, "too many curves.")
        self.failUnless(l01 in cc, "l01 is missing")
        self.failUnless(l12 in cc, "l12 is missing")
        self.failUnless(l20 in cc, "l20 is missing")

        p=c.getPrimitives()
        self.failUnless(len(p) == 9, "too many primitives.")
        self.failUnless(l01 in p, "l01 is missing")
        self.failUnless(l12 in p, "l21 is missing")
        self.failUnless(l20 in p, "l20 is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")

        cp=c.copy()
        self.failUnless(isinstance(cp,CurveLoop), "copy returns is not an arc.")
        self.failUnless(not cp == c, "copy equals source")
        self.failUnless(cp.isColocated(c),"copy is not colocated with its source.")
        cc=cp.getCurves()
        self.failUnless(len(cc) == 3, "too many primitives in copy.")
        self.failUnless(not l01 in cc,"copy uses l01.")
        self.failUnless(not l12 in cc,"copy uses l12.")
        self.failUnless(not l20 in cc,"copy uses l20.")
         
        p0_m=Point(0,0,0)
        p1_m=Point(-1,-1,-1)
        p2_m=Point(-2,-2,-2)
        p3_m=Point(-3,-3,-3)
        p4_m=Point(-1,-2,-3)

        l01_m=Line(p0_m,p1_m)
        l12_m=Arc(p3_m,p1_m,p2_m)
        l20_m=Spline(p2_m,p4_m,p0_m)

        dc=c.apply(Dilation(-1.))
        self.failUnless(dc.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"dilation is wrong.")
        cc=dc.getCurves()
        self.failUnless(len(cc) == 3, "too many primitives in dilation result.")
        self.failUnless(not l01 in cc,"l01 is in dilation result.")
        self.failUnless(not l12 in cc,"l12 is in dilation result.")
        self.failUnless(not l20 in cc,"l20 is in dilation result.")

        c.modifyBy(Dilation(-1.))
        self.failUnless(c.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"inplace dilation is wrong.")
        cc=c.getCurves()
        self.failUnless(len(cc) == 3, "too many primitives in modified object.")
        self.failUnless(l01 in cc,"l01 missed in  modified object.")
        self.failUnless(cc[cc.index(l01)].hasSameOrientation(l01),"l01 in modified object has wrong orientation.")
        self.failUnless(l12 in cc,"l12 missed in  modified object.")
        self.failUnless(cc[cc.index(l12)].hasSameOrientation(l12),"l12 in modified object has wrong orientation.")
        self.failUnless(l20 in cc,"l20 missed in  modified object.")
        self.failUnless(cc[cc.index(l20)].hasSameOrientation(l20),"l20 in modified object has wrong orientation.")
      
   def test_ReverseCurveLoop(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        p5=Point(10,20,3)
        p6=Point(1,2,30)

        l01=Line(p0,p1)
        l12=Arc(p3,p1,p2)
        l20=Spline(p2,p4,p0)

        lx=Line(p2,p3)
        ly=Line(p3,p1)

        CC0=CurveLoop(l01,l20,l12)
        c=-CC0

        self.failUnless(c.hasSameOrientation(c),"has not same orientation like itself")
        self.failUnless(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.failUnless(not c.isColocated(p4),"-CurveLoop is colocated with point.")
        self.failUnless(c.isColocated(c),"-CurveLoop is not colocated with its self.")
        self.failUnless(c.isColocated(CurveLoop(l01,l12,l20)),"-CurveLoop is not colocated with its copy.")
        self.failUnless(c.isColocated(CurveLoop(l20,l01,l12)),"-CurveLoop is not colocated with its copy with shifted points.")
        self.failUnless(c.isColocated(CurveLoop(l20,l12,l01)),"-CurveLoop is not colocated with its copy with shuffled points.")
        self.failUnless(not c.isColocated(CurveLoop(lx,ly,l12)),"-CurveLoop is colocated with different CurveLoop.")

        self.failUnless(len(c) == 3, "wrong length")

        c.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")


        cc=c.getCurves()
        self.failUnless(len(cc) == 3, "too many curves.")
        self.failUnless(l01 in cc, "l01 is missing")
        self.failUnless(l12 in cc, "l12 is missing")
        self.failUnless(l20 in cc, "l20 is missing")

        p=c.getPrimitives()
        self.failUnless(len(p) == 9, "too many primitives.")
        self.failUnless(l01 in p, "l01 is missing")
        self.failUnless(l12 in p, "l21 is missing")
        self.failUnless(l20 in p, "l20 is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")

        cp=c.copy()
        self.failUnless(isinstance(cp,ReverseCurveLoop), "copy returns is not an ReverseCurveLoop.")
        self.failUnless(not cp == c, "copy equals source")
        self.failUnless(cp.isColocated(c),"copy is not colocated with its source.")
        cc=cp.getCurves()
        self.failUnless(len(cc) == 3, "too many primitives in copy.")
        self.failUnless(not l01 in cc,"copy uses l01.")
        self.failUnless(not l12 in cc,"copy uses l12.")
        self.failUnless(not l20 in cc,"copy uses l20.")
         
        p0_m=Point(0,0,0)
        p1_m=Point(-1,-1,-1)
        p2_m=Point(-2,-2,-2)
        p3_m=Point(-3,-3,-3)
        p4_m=Point(-1,-2,-3)

        l01_m=Line(p0_m,p1_m)
        l12_m=Arc(p3_m,p1_m,p2_m)
        l20_m=Spline(p2_m,p4_m,p0_m)

        dc=c.apply(Dilation(-1.))
        self.failUnless(dc.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"dilation is wrong.")
        cc=dc.getCurves()
        self.failUnless(len(cc) == 3, "too many primitives in dilation result.")
        self.failUnless(not l01 in cc,"l01 is in dilation result.")
        self.failUnless(not l12 in cc,"l12 is in dilation result.")
        self.failUnless(not l20 in cc,"l20 is in dilation result.")

        c.modifyBy(Dilation(-1.))
        self.failUnless(c.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"inplace dilation is wrong.")
        cc=c.getCurves()
        self.failUnless(len(cc) == 3, "too many primitives in modified object.")
        self.failUnless(l01 in cc,"l01 missed in  modified object.")
        self.failUnless(cc[cc.index(l01)].hasSameOrientation(-l01),"l01 in modified object has wrong orientation.")
        self.failUnless(l12 in cc,"l12 missed in  modified object.")
        self.failUnless(cc[cc.index(l12)].hasSameOrientation(-l12),"l12 in modified object has wrong orientation.")
        self.failUnless(l20 in cc,"l20 missed in  modified object.")
        self.failUnless(cc[cc.index(l20)].hasSameOrientation(-l20),"l20 in modified object has wrong orientation.")
      
   def test_RuledSurface(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        p5=Point(10,20,3)
        p6=Point(1,2,30)

        l01=Line(p0,p1)
        l12_1=Arc(p3,p1,p2)
        l12_2_1=Spline(p1,p3,p4)
        l12_2_2=Spline(p4,p5,p2)
        l12_3=Line(p1,p2)
        l20=Spline(p2,p4,p0)

        cl1=CurveLoop(l01,l12_1,l20) 
        cl2=CurveLoop(l01,l12_2_1,l12_2_2,l20)
        cl3=CurveLoop(l01,l12_3,l20)

        self.failUnlessRaises(TypeError,RuledSurface,l01)

        s=RuledSurface(cl1)

        cl=s.getBoundaryLoop()
        self.failUnless(cl == cl1, " wrong boundary loops")
        self.failUnless(cl.hasSameOrientation(cl1),"cl1 has incorrect orientation.")

        self.failUnless(s.hasSameOrientation(s),"has not same orientation like itself")
        self.failUnless(not s.hasSameOrientation(-s),"has same orientation like -itself")

        crvs=s.getBoundary()
        self.failUnless(len(crvs) == 3, "too many boundary corves.")
        self.failUnless(l01 in crvs, "l01 is missing in boundary")
        self.failUnless(crvs[crvs.index(l01)].hasSameOrientation(l01),"l01 has incorrect orientation.")
        self.failUnless(l12_1 in crvs, "l21 is missing in boundary")
        self.failUnless(crvs[crvs.index(l12_1)].hasSameOrientation(l12_1),"l12_1 has incorrect orientation.")
        self.failUnless(l20 in crvs, "l20 is missing in boundary")
        self.failUnless(crvs[crvs.index(l20)].hasSameOrientation(l20),"l12_1 has incorrect orientation.")
               

        self.failUnless(not s.isColocated(p4),"RuledSurface is colocated with point.")
        self.failUnless(s.isColocated(s),"RuledSurface is not colocated with its self.")
        self.failUnless(s.isColocated(RuledSurface(cl1)),"RuledSurface is not colocated with its copy.")
        self.failUnless(not s.isColocated(RuledSurface(cl2)),"RuledSurface is colocated with different length")
        self.failUnless(not s.isColocated(RuledSurface(cl3)),"RuledSurface is colocated with same length.")

        s.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")

        p=s.getPrimitives()
        self.failUnless(len(p) == 10, "too many primitives.")
        self.failUnless(cl1 in p, "cl1 is missing")
        self.failUnless(l01 in p, "l01 is missing")
        self.failUnless(l12_1 in p, "l21 is missing")
        self.failUnless(l20 in p, "l20 is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")

        sp=s.copy()
        self.failUnless(isinstance(sp,RuledSurface), "copy returns is not a RuledSurface.")
        self.failUnless(not sp == s, "copy equals source")
        self.failUnless(sp.isColocated(s),"copy is not colocated with its source.")
        cbl=sp.getBoundaryLoop()
        self.failUnless(not cbl == cl1,"copy uses cl1.")
        cp=sp.getPrimitives()
        self.failUnless(len(cp) == 10, "copy as too many primitives.")
        self.failUnless(not cl1 in cp, "copy is using cl1")
        self.failUnless(not l01 in cp, "copy is using l01")
        self.failUnless(not l12_1 in cp, "copy is using l21")
        self.failUnless(not l20 in cp, "copy is using l20")
        self.failUnless(not p0 in cp, "copy is using p0")
        self.failUnless(not p1 in cp, "copy is using p1")
        self.failUnless(not p2 in cp, "copy is using p2")
        self.failUnless(not p3 in cp, "copy is using p3")
        self.failUnless(not p4 in cp, "copy is using p4")
        del cp
         
        p0_m=Point(0,0,0)
        p1_m=Point(-1,-1,-1)
        p2_m=Point(-2,-2,-2)
        p3_m=Point(-3,-3,-3)
        p4_m=Point(-1,-2,-3)

        l01_m=Line(p0_m,p1_m)
        l12_m=Arc(p3_m,p1_m,p2_m)
        l20_m=Spline(p2_m,p4_m,p0_m)

        ds=s.apply(Dilation(-1.))
        self.failUnless(ds.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"dilation is wrong.")
        cbl=ds.getBoundaryLoop()
        self.failUnless(not cbl == cl1,"dilation uses cl1.")
        cp=ds.getPrimitives()
        self.failUnless(len(cp) == 10, "dilation as too many primitives.")
        self.failUnless(not cl1 in cp, "dilation is using cl1")
        self.failUnless(not l01 in cp, "dilation is using l01")
        self.failUnless(not l12_1 in cp, "dilation is using l21")
        self.failUnless(not l20 in cp, "dilation is using l20")
        self.failUnless(not p0 in cp, "dilation is using p0")
        self.failUnless(not p1 in cp, "dilation is using p1")
        self.failUnless(not p2 in cp, "dilation is using p2")
        self.failUnless(not p3 in cp, "dilation is using p3")
        self.failUnless(not p4 in cp, "dilation is using p4")

        s.modifyBy(Dilation(-1.))
        self.failUnless(s.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"inplace dilation is wrong.")

        p=s.getPrimitives()
        self.failUnless(len(p) == 10, "inplace dilation has too many primitives.")
        self.failUnless(cl1 in p, "inplace dilation cl1 is missing")
        self.failUnless(l01 in p, "inplace dilation l01 is missing")
        self.failUnless(l12_1 in p, "inplace dilation l21 is missing")
        self.failUnless(l20 in p, "inplace dilation l20 is missing")
        self.failUnless(p0 in p, "inplace dilation p0 is missing")
        self.failUnless(p1 in p, "inplace dilation p1 is missing")
        self.failUnless(p2 in p, "inplace dilation p2 is missing")
        self.failUnless(p3 in p, "inplace dilation p3 is missing")
        self.failUnless(p4 in p, "inplace dilation p4 is missing")

        p=s.getBoundary()
        self.failUnless(len(p) == 3, "inplace dilation has too many boundary curves.")
        self.failUnless(l01 in p, "inplace dilation l01 is missing in boundary curves.")
        self.failUnless(p[p.index(l01)].hasSameOrientation(l01),"l01 in getBoundary after dilation has incorrect orientation.")
        self.failUnless(l12_1 in p, "inplace dilation l21 is missing")
        self.failUnless(p[p.index(l12_1)].hasSameOrientation(l12_1),"l12_1 in getBoundary after dilation has incorrect orientation.")
        self.failUnless(l20 in p, "inplace dilation l20 is missing")
        self.failUnless(p[p.index(l20)].hasSameOrientation(l20),"l20 in getBoundary after dilation has incorrect orientation.")

        p=s.getBoundaryLoop()
        self.failUnless(cl1 == p, "inplace dilation s.getBoundaryLoop does not return cl1")
        self.failUnless(p.hasSameOrientation(cl1),"cl1 in getBoundaryLoop after dilation has incorrect orientation.")

   def test_ReverseRuledSurface(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        p5=Point(10,20,3)
        p6=Point(1,2,30)

        l01=Line(p0,p1)
        l12_1=Arc(p3,p1,p2)
        l12_2_1=Spline(p1,p3,p4)
        l12_2_2=Spline(p4,p5,p2)
        l12_3=Line(p1,p2)
        l20=Spline(p2,p4,p0)

        cl1=CurveLoop(l01,l12_1,l20) 
        cl2=CurveLoop(l01,l12_2_1,l12_2_2,l20)
        cl3=CurveLoop(l01,l12_3,l20)

        self.failUnlessRaises(TypeError,RuledSurface,l01)

        CC0=RuledSurface(cl1)
        s=-CC0

        cl=s.getBoundaryLoop()
        self.failUnless(cl == cl1, " wrong boundary loops")
        self.failUnless(cl.hasSameOrientation(-cl1),"cl1 has incorrect orientation.")

        self.failUnless(s.hasSameOrientation(s),"has not same orientation like itself")
        self.failUnless(not s.hasSameOrientation(-s),"has same orientation like -itself")

        crvs=s.getBoundary()
        self.failUnless(len(crvs) == 3, "too many boundary corves.")
        self.failUnless(l01 in crvs, "l01 is missing in boundary")
        self.failUnless(crvs[crvs.index(l01)].hasSameOrientation(-l01),"l01 has incorrect orientation.")
        self.failUnless(l12_1 in crvs, "l21 is missing in boundary")
        self.failUnless(crvs[crvs.index(l12_1)].hasSameOrientation(-l12_1),"l12_1 has incorrect orientation.")
        self.failUnless(l20 in crvs, "l20 is missing in boundary")
        self.failUnless(crvs[crvs.index(l20)].hasSameOrientation(-l20),"l12_1 has incorrect orientation.")
               

        self.failUnless(not s.isColocated(p4),"RuledSurface is colocated with point.")
        self.failUnless(s.isColocated(s),"RuledSurface is not colocated with its self.")
        self.failUnless(s.isColocated(RuledSurface(cl1)),"RuledSurface is not colocated with its copy.")
        self.failUnless(not s.isColocated(RuledSurface(cl2)),"RuledSurface is colocated with different length")
        self.failUnless(not s.isColocated(RuledSurface(cl3)),"RuledSurface is colocated with same length.")

        s.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")

        p=s.getPrimitives()
        self.failUnless(len(p) == 10, "too many primitives.")
        self.failUnless(cl1 in p, "cl1 is missing")
        self.failUnless(l01 in p, "l01 is missing")
        self.failUnless(l12_1 in p, "l21 is missing")
        self.failUnless(l20 in p, "l20 is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")

        sp=s.copy()
        self.failUnless(isinstance(sp,ReverseRuledSurface), "copy returns is not a RuledSurface.")
        self.failUnless(not sp == s, "copy equals source")
        self.failUnless(sp.isColocated(s),"copy is not colocated with its source.")
        cbl=sp.getBoundaryLoop()
        self.failUnless(not cbl == cl1,"copy uses cl1.")
        cp=sp.getPrimitives()
        self.failUnless(len(cp) == 10, "copy as too many primitives.")
        self.failUnless(not cl1 in cp, "copy is using cl1")
        self.failUnless(not l01 in cp, "copy is using l01")
        self.failUnless(not l12_1 in cp, "copy is using l21")
        self.failUnless(not l20 in cp, "copy is using l20")
        self.failUnless(not p0 in cp, "copy is using p0")
        self.failUnless(not p1 in cp, "copy is using p1")
        self.failUnless(not p2 in cp, "copy is using p2")
        self.failUnless(not p3 in cp, "copy is using p3")
        self.failUnless(not p4 in cp, "copy is using p4")
        del cp
         
        p0_m=Point(0,0,0)
        p1_m=Point(-1,-1,-1)
        p2_m=Point(-2,-2,-2)
        p3_m=Point(-3,-3,-3)
        p4_m=Point(-1,-2,-3)

        l01_m=Line(p0_m,p1_m)
        l12_m=Arc(p3_m,p1_m,p2_m)
        l20_m=Spline(p2_m,p4_m,p0_m)

        ds=s.apply(Dilation(-1.))
        self.failUnless(ds.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"dilation is wrong.")
        cbl=ds.getBoundaryLoop()
        self.failUnless(not cbl == cl1,"dilation uses cl1.")
        cp=ds.getPrimitives()
        self.failUnless(len(cp) == 10, "dilation as too many primitives.")
        self.failUnless(not cl1 in cp, "dilation is using cl1")
        self.failUnless(not l01 in cp, "dilation is using l01")
        self.failUnless(not l12_1 in cp, "dilation is using l21")
        self.failUnless(not l20 in cp, "dilation is using l20")
        self.failUnless(not p0 in cp, "dilation is using p0")
        self.failUnless(not p1 in cp, "dilation is using p1")
        self.failUnless(not p2 in cp, "dilation is using p2")
        self.failUnless(not p3 in cp, "dilation is using p3")
        self.failUnless(not p4 in cp, "dilation is using p4")

        s.modifyBy(Dilation(-1.))
        self.failUnless(s.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"inplace dilation is wrong.")

        p=s.getPrimitives()
        self.failUnless(len(p) == 10, "inplace dilation has too many primitives.")
        self.failUnless(cl1 in p, "inplace dilation cl1 is missing")
        self.failUnless(l01 in p, "inplace dilation l01 is missing")
        self.failUnless(l12_1 in p, "inplace dilation l21 is missing")
        self.failUnless(l20 in p, "inplace dilation l20 is missing")
        self.failUnless(p0 in p, "inplace dilation p0 is missing")
        self.failUnless(p1 in p, "inplace dilation p1 is missing")
        self.failUnless(p2 in p, "inplace dilation p2 is missing")
        self.failUnless(p3 in p, "inplace dilation p3 is missing")
        self.failUnless(p4 in p, "inplace dilation p4 is missing")

        p=s.getBoundary()
        self.failUnless(len(p) == 3, "inplace dilation has too many boundary curves.")
        self.failUnless(l01 in p, "inplace dilation l01 is missing in boundary curves.")
        self.failUnless(p[p.index(l01)].hasSameOrientation(-l01),"l01 in getBoundary after dilation has incorrect orientation.")
        self.failUnless(l12_1 in p, "inplace dilation l21 is missing")
        self.failUnless(p[p.index(l12_1)].hasSameOrientation(-l12_1),"l12_1 in getBoundary after dilation has incorrect orientation.")
        self.failUnless(l20 in p, "inplace dilation l20 is missing")
        self.failUnless(p[p.index(l20)].hasSameOrientation(-l20),"l20 in getBoundary after dilation has incorrect orientation.")

        p=s.getBoundaryLoop()
        self.failUnless(cl1 == p, "inplace dilation s.getBoundaryLoop does not return cl1")
        self.failUnless(p.hasSameOrientation(-cl1),"cl1 in getBoundaryLoop after dilation has incorrect orientation.")


   def test_PlaneSurface(self):
        p0=Point(0,0,0,0.1)
        p1=Point(10,0,0,0.2)
        p2=Point(10,10,0,0.3)
        p3=Point(0,10,3,0.4)
        p4=Point(5,5,0,0.001)
        p5=Point(7,5,0,0.001)
        p6=Point(5,7,0,0.001)
        p7=Point(8,8,0,0.001)
        p8=Point(9,9,0,0.001)

        l0=Line(p0,p1)
        l1=Line(p1,p2)
        l2=Line(p2,p3)
        l3=Line(p3,p0)

        l9=Line(p1,p8)
        l10=Line(p8,p3)

        l4=Line(p4,p5)
        l5=Line(p5,p6)
        l6=Line(p6,p4)
        l7=Line(p6,p7)
        l8=Line(p7,p4)

        a1=Arc(p4,p3,p1)
        a2=Arc(p7,p5,p6)

        cl=CurveLoop(l0,l1,l2,l3)
        h=CurveLoop(l4,l5,l6)
        cl_s=CurveLoop(l0,l9,l10,l3)
        h2=CurveLoop(l4,l5,l7,l8)
        cl_a=CurveLoop(a1,l1,l2)
        h_a=CurveLoop(a2,l6,l4)

        self.failUnlessRaises(TypeError,PlaneSurface,l4)
        self.failUnlessRaises(TypeError,PlaneSurface,cl_a,h)
        self.failUnlessRaises(TypeError,PlaneSurface,cl,[h_a])

        s=PlaneSurface(cl,holes=[h])

        cl2=s.getBoundaryLoop()
        self.failUnless(cl == cl2, " wrong boundary loops")
        self.failUnless(cl.hasSameOrientation(cl2),"cl has incorrect orientation.")

        hs=s.getHoles()
        self.failUnless(len(hs) == 1, "one holes expected.")
        self.failUnless(h==hs[0], "h is not defined as hole.")
        self.failUnless(hs[0].hasSameOrientation(h),"hs has incorrect orientation.")

        self.failUnless(s.hasSameOrientation(s),"has not same orientation like itself")
        self.failUnless(not s.hasSameOrientation(-s),"has same orientation like -itself")

        crvs=s.getBoundary()
        self.failUnless(len(crvs) == 7, "too many boundary corves.")
        self.failUnless(l0 in crvs, "l0 is missing in boundary")
        self.failUnless(crvs[crvs.index(l0)].hasSameOrientation(l0),"l0 has incorrect orientation.")
        self.failUnless(l1 in crvs, "l1 is missing in boundary")
        self.failUnless(crvs[crvs.index(l1)].hasSameOrientation(l1),"l1 has incorrect orientation.")
        self.failUnless(l2 in crvs, "l2 is missing in boundary")
        self.failUnless(crvs[crvs.index(l2)].hasSameOrientation(l2),"l2 has incorrect orientation.")
        self.failUnless(l3 in crvs, "l3 is missing in boundary")
        self.failUnless(crvs[crvs.index(l3)].hasSameOrientation(l3),"l3 has incorrect orientation.")
        self.failUnless(l4 in crvs, "l4 is missing in boundary")
        self.failUnless(crvs[crvs.index(l4)].hasSameOrientation(l4),"l4 has incorrect orientation.")
        self.failUnless(l5 in crvs, "l5 is missing in boundary")
        self.failUnless(crvs[crvs.index(l5)].hasSameOrientation(l5),"l5 has incorrect orientation.")
        self.failUnless(l6 in crvs, "l6 is missing in boundary")
        self.failUnless(crvs[crvs.index(l6)].hasSameOrientation(l6),"l6 has incorrect orientation.")
               
        self.failUnless(not s.isColocated(p4),"PlaneSurface is colocated with point.")
        self.failUnless(s.isColocated(s),"PlaneSurface is not colocated with its self.")
        self.failUnless(s.isColocated(PlaneSurface(cl,holes=[h])),"PlaneSurface is not colocated with its copy.")
        self.failUnless(not s.isColocated(PlaneSurface(cl)),"PlaneSurface is colocated with PlaneSurface with same boundary but no hole")
        self.failUnless(not s.isColocated(PlaneSurface(cl_s,holes=[h])),"PlaneSurface is colocated with PlaneSurface with deformed boundary")
        self.failUnless(not s.isColocated(PlaneSurface(cl,holes=[h2])),"PlaneSurface is colocated with modified hole")

        s.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.failUnless(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.failUnless(p6.getLocalScale()==3., "p6 has wrong local scale.")

        p=s.getPrimitives()
        self.failUnless(len(p) == 17, "too many primitives.")
        self.failUnless(s in p, "cl is missing")
        self.failUnless(cl in p, "cl is missing")
        self.failUnless(h in p, "h is missing")
        self.failUnless(l0 in p, "l0 is missing")
        self.failUnless(l1 in p, "l1 is missing")
        self.failUnless(l2 in p, "l2 is missing")
        self.failUnless(l3 in p, "l3 is missing")
        self.failUnless(l4 in p, "l4 is missing")
        self.failUnless(l5 in p, "l5 is missing")
        self.failUnless(l6 in p, "l6 is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")
        self.failUnless(p5 in p, "p5 is missing")
        self.failUnless(p6 in p, "p6 is missing")

        sp=s.copy()
        self.failUnless(isinstance(sp,PlaneSurface), "copy returns is not a PlaneSurface.")
        self.failUnless(not sp == s, "copy equals source")
        self.failUnless(sp.isColocated(s),"copy is not colocated with its source.")
        cbl=sp.getBoundaryLoop()
        self.failUnless(not cbl == cl,"copy uses cl1.")
        hs=sp.getHoles()
        self.failUnless(len(hs)==1,"copy is missing holes.")
        self.failUnless(not hs[0]== h,"copy contains h as a hole.")
        cp=sp.getPrimitives()
        self.failUnless(len(cp) == 17, "copy as too many primitives.")
        self.failUnless(not s in cp, "copy contains s")
        self.failUnless(not cl in cp, "copy contains cl")
        self.failUnless(not h in cp, "copy contains h")
        self.failUnless(not l0 in cp, "copy contains l0")
        self.failUnless(not l1 in cp, "copy contains l1")
        self.failUnless(not l2 in cp, "copy contains l2")
        self.failUnless(not l3 in cp, "copy contains l3")
        self.failUnless(not l4 in cp, "copy contains l4")
        self.failUnless(not l5 in cp, "copy contains l5")
        self.failUnless(not l6 in cp, "copy contains l6")
        self.failUnless(not p0 in cp, "copy contains p0")
        self.failUnless(not p1 in cp, "copy contains p1")
        self.failUnless(not p2 in cp, "copy contains p2")
        self.failUnless(not p3 in cp, "copy contains p3")
        self.failUnless(not p4 in cp, "copy contains p4")
        self.failUnless(not p5 in cp, "copy contains p5")
        self.failUnless(not p6 in cp, "copy contains p6")
        del sp
         

        p0_m=Point(0,0,0,0.1)
        p1_m=Point(-10,0,0,0.2)
        p2_m=Point(-10,-10,0,0.3)
        p3_m=Point(0,-10,-3,0.4)
        p4_m=Point(-5,-5,0,0.001)
        p5_m=Point(-7,-5,0,0.001)
        p6_m=Point(-5,-7,0,0.001)

        l0_m=Line(p0_m,p1_m)
        l1_m=Line(p1_m,p2_m)
        l2_m=Line(p2_m,p3_m)
        l3_m=Line(p3_m,p0_m)

        l4_m=Line(p4_m,p5_m)
        l5_m=Line(p5_m,p6_m)
        l6_m=Line(p6_m,p4_m)

        cl_m=CurveLoop(l0_m,l1_m,l2_m,l3_m)
        h_m=CurveLoop(l4_m,l5_m,l6_m)

        ds=s.apply(Dilation(-1.))
        self.failUnless(ds.isColocated(PlaneSurface(cl_m,holes=[h_m])),"dilation is wrong.")
        cbl=ds.getBoundaryLoop()
        self.failUnless(not cbl == cl,"dilation uses cl1.")
        hs=ds.getHoles()
        self.failUnless(len(hs)==1,"dilation is missing holes.")
        self.failUnless(not hs[0]== h,"dilation contains h as a hole.")
        cp=ds.getPrimitives()
        self.failUnless(len(cp) == 17, "dilation as too many primitives.")
        self.failUnless(not s in cp, "dilation contains s")
        self.failUnless(not cl in cp, "dilation contains cl")
        self.failUnless(not h in cp, "dilation contains h")
        self.failUnless(not l0 in cp, "dilation contains l0")
        self.failUnless(not l1 in cp, "dilation contains l1")
        self.failUnless(not l2 in cp, "dilation contains l2")
        self.failUnless(not l3 in cp, "dilation contains l3")
        self.failUnless(not l4 in cp, "dilation contains l4")
        self.failUnless(not l5 in cp, "dilation contains l5")
        self.failUnless(not l6 in cp, "dilation contains l6")
        self.failUnless(not p0 in cp, "dilation contains p0")
        self.failUnless(not p1 in cp, "dilation contains p1")
        self.failUnless(not p2 in cp, "dilation contains p2")
        self.failUnless(not p3 in cp, "dilation contains p3")
        self.failUnless(not p4 in cp, "dilation contains p4")
        self.failUnless(not p5 in cp, "dilation contains p5")
        self.failUnless(not p6 in cp, "dilation contains p6")

        s.modifyBy(Dilation(-1.))
        self.failUnless(s.isColocated(PlaneSurface(cl_m,holes=[h_m])),"inplace dilation is wrong.")
        cbl=s.getBoundaryLoop()
        self.failUnless(cbl == cl,"inplace dilation does not use cl1.")
        self.failUnless(cl.hasSameOrientation(cbl),"cl has incorrect orientation.")
        hs=s.getHoles()
        self.failUnless(len(hs)==1,"inplace dilation is missing holes.")
        self.failUnless(hs[0]== h,"inplace dilation must contain h as a hole.")
        self.failUnless(hs[0].hasSameOrientation(h),"hole in inplace dilation has incorrect orientation.")

        cp=s.getPrimitives()
        self.failUnless(len(cp) == 17, "inplace dilation as too many primitives.")
        self.failUnless(s in cp, "inplace dilation must use s")
        self.failUnless(cl in cp, "inplace dilation must use cl")
        self.failUnless(h in cp, "inplace dilation must use h")
        self.failUnless(l0 in cp, "inplace dilation must use l0")
        self.failUnless(l1 in cp, "inplace dilation must use l1")
        self.failUnless(l2 in cp, "inplace dilation must use l2")
        self.failUnless(l3 in cp, "inplace dilation must use l3")
        self.failUnless(l4 in cp, "inplace dilation must use l4")
        self.failUnless(l5 in cp, "inplace dilation must use l5")
        self.failUnless(l6 in cp, "inplace dilation must use l6")
        self.failUnless(p0 in cp, "inplace dilation must use p0")
        self.failUnless(p1 in cp, "inplace dilation must use p1")
        self.failUnless(p2 in cp, "inplace dilation must use p2")
        self.failUnless(p3 in cp, "inplace dilation must use p3")
        self.failUnless(p4 in cp, "inplace dilation must use p4")
        self.failUnless(p5 in cp, "inplace dilation must use p5")
        self.failUnless(p6 in cp, "inplace dilation must use p6")

   def test_SurfaceLoop(self):
        p0=Point( 0, 0, 0,0.1)
        p1=Point(10, 0, 0,0.1)
        p2=Point( 0,10, 0,0.1)
        p3=Point(10,10, 0,0.1)
        p4=Point( 0, 0,10,0.1)
        p5=Point(10, 0,10,0.1)
        p6=Point( 0,10,10,0.1)
        p7=Point(10,10,10,0.1)

        q0=Point( 4, 0, 4,0.1)
        q1=Point( 6, 0, 4,0.1)
        q2=Point( 4,10, 4,0.1)
        q3=Point( 6,10, 4,0.1)
        q4=Point( 4, 0, 6,0.1)
        q5=Point( 6, 0, 6,0.1)
        q6=Point( 4,10, 6,0.1)
        q7=Point( 6,10, 6,0.1)

        l01=Line(p0,p1)
        l15=Line(p1,p5)
        l54=Line(p5,p4)
        l04=Line(p4,p0)

        l13=Line(p1,p3)
        l37=Line(p3,p7)
        l75=Line(p7,p5)
        l67=Line(p6,p7)
        l26=Line(p2,p6)
        l32=Line(p3,p2)
        l20=Line(p2,p0)
        l46=Line(p4,p6)

        m01=Line(q0,q1)
        m15=Line(q1,q5)
        m54=Line(q5,q4)
        m40=Line(q4,q0)
        m23=Line(q2,q3)
        m37=Line(q3,q7)
        m76=Line(q7,q6)
        m62=Line(q6,q2)

        m02=Line(q0,q2)
        m13=Line(q1,q3)
        m46=Line(q4,q6)
        m57=Line(q5,q7)

        cl_l1=CurveLoop(l01,l15,l54,l04)
        cl_m1=CurveLoop(m01,m15,m54,m40)
        s1=PlaneSurface(cl_l1,holes=[cl_m1])
        s1_v=PlaneSurface(cl_l1)

        cl_l2=CurveLoop(-l15,l13,l37,l75)
        s2=PlaneSurface(cl_l2)

        cl_l3=CurveLoop(l32,-l37,l67,l26)
        cl_m3=CurveLoop(-m23,-m37,-m76,-m62)
        s3=PlaneSurface(cl_l3,holes=[cl_m3])
        s3_v=PlaneSurface(cl_l3)
     
     
        cl_l4=CurveLoop(l20,-l26,l46,-l04)
        s4=PlaneSurface(cl_l4)

        cl_l5=CurveLoop(l32,l20,l01,l13)
        s5=PlaneSurface(-cl_l5)

        cl_l6=CurveLoop(l67,l75,l54,l46)
        s6=PlaneSurface(-cl_l6)

        cl_m7=CurveLoop(m13,m37,-m57,-m15)
        s7=PlaneSurface(cl_m7)
        
        cl_m8=CurveLoop(m57,m76,-m46,-m54)
        s8=PlaneSurface(cl_m8)

        cl_m9=CurveLoop(m46,m62,-m02,-m40)
        s9=PlaneSurface(cl_m9)

        cl_m10=CurveLoop(-m01,m02,m23,-m13)
        s10=PlaneSurface(cl_m10)

        # self.failUnlessRaises(ValueError,SurfaceLoop,s1,s3)
        # self.failUnlessRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5)
        # self.failUnlessRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5,s5)
        s=SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)

        self.failUnless(s.hasSameOrientation(s),"has not same orientation like itself")
        self.failUnless(not s.hasSameOrientation(-s),"has same orientation like -itself")

        cc=s.getSurfaces()
        self.failUnless(len(cc) == 10, "too many curves.")
        self.failUnless(s1 in cc, "s1 is missing")
        self.failUnless(s2 in cc, "s2 is missing")
        self.failUnless(s3 in cc, "s3 is missing")
        self.failUnless(s4 in cc, "s4 is missing")
        self.failUnless(s5 in cc, "s5 is missing")
        self.failUnless(s6 in cc, "s6 is missing")
        self.failUnless(s7 in cc, "s7 is missing")
        self.failUnless(s8 in cc, "s8 is missing")
        self.failUnless(s9 in cc, "s9 is missing")
        self.failUnless(s10 in cc, "s10 is missing")

        self.failUnless(not s.isColocated(p4),"SurfaceLoop is colocated with point.")
        self.failUnless(s.isColocated(s),"SurfaceLoop is not colocated with its self.")
        self.failUnless(s.isColocated(-s),"SurfaceLoop is not colocated with its reverse.")
        self.failUnless(s.isColocated(SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)),"SurfaceLoop is not colocated with its copy.")
        self.failUnless(s.isColocated(SurfaceLoop(-s10, s1,s2,s3,s4,s5,s6,-s7,-s8,-s9)),"SurfaceLoop is not colocated with its copy with shifted points.")
        self.failUnless(s.isColocated(SurfaceLoop(s1, -s7, s2, -s9, s3, -s8, s4, -s10, s5,s6)),"SurfaceLoop is not colocated with its copy with shuffled points.")
        self.failUnless(not s.isColocated(SurfaceLoop(s1_v,s2,s3_v,s4,s5,s6)),"SurfaceLoop is colocated with different SurfaceLoop.")

        self.failUnless(len(s) == 10, "wrong length")

        s.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.failUnless(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.failUnless(p6.getLocalScale()==3., "p6 has wrong local scale.")
        self.failUnless(p7.getLocalScale()==3., "p7 has wrong local scale.")

        p=s.getPrimitives()
        self.failUnless(len(p) == 63, "too many primitives.")
        self.failUnless(s in p, "s is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")
        self.failUnless(p5 in p, "p5 is missing")
        self.failUnless(p6 in p, "p6 is missing")
        self.failUnless(p7 in p, "p7 is missing")
        self.failUnless(q0 in p, "q0 is missing")
        self.failUnless(q1 in p, "q1 is missing")
        self.failUnless(q2 in p, "q2 is missing")
        self.failUnless(q3 in p, "q3 is missing")
        self.failUnless(q4 in p, "q4 is missing")
        self.failUnless(q5 in p, "q5 is missing")
        self.failUnless(q6 in p, "q6 is missing")
        self.failUnless(q7 in p, "q7 is missing")
        self.failUnless(l01 in p, "l01 is missing")
        self.failUnless(l15 in p, "l15 is missing")
        self.failUnless(l54 in p, "l54 is missing")
        self.failUnless(l04 in p, "l04 is missing")
        self.failUnless(l13 in p, "l13 is missing")
        self.failUnless(l37 in p, "l37 is missing")
        self.failUnless(l75 in p, "l75 is missing")
        self.failUnless(l67 in p, "l67 is missing")
        self.failUnless(l26 in p, "l26 is missing")
        self.failUnless(l32 in p, "l32 is missing")
        self.failUnless(l20 in p, "l20 is missing")
        self.failUnless(l46 in p, "l46 is missing")
        self.failUnless(m01 in p, "m01 is missing")
        self.failUnless(m15 in p, "m15 is missing")
        self.failUnless(m54 in p, "m54 is missing")
        self.failUnless(m40 in p, "m40 is missing")
        self.failUnless(m23 in p, "m23 is missing")
        self.failUnless(m37 in p, "m37 is missing")
        self.failUnless(m76 in p, "m76 is missing")
        self.failUnless(m62 in p, "m62 is missing")
        self.failUnless(m02 in p, "m02 is missing")
        self.failUnless(m13 in p, "m13 is missing")
        self.failUnless(m46 in p, "m46 is missing")
        self.failUnless(m57 in p, "m57 is missing")
        self.failUnless(cl_l1 in p, "cl_l1 is missing")
        self.failUnless(cl_m1 in p, "cl_m1 is missing")
        self.failUnless(s1 in p, "s1 is missing")
        self.failUnless(cl_l2 in p, "cl_l2 is missing")
        self.failUnless(s2 in p, "s2 is missing")
        self.failUnless(cl_l3 in p, "cl_l3 is missing")
        self.failUnless(cl_m3 in p, "cl_m3 is missing")
        self.failUnless(s3 in p, "s3 is missing")
        self.failUnless(cl_l4 in p, "cl_l4 is missing")
        self.failUnless(s4 in p, "s4 is missing")
        self.failUnless(cl_l5 in p, "cl_l5 is missing")
        self.failUnless(s5 in p, "s5 is missing")
        self.failUnless(cl_l6 in p, "cl_l6 is missing")
        self.failUnless(s6 in p, "s6 is missing")
        self.failUnless(cl_m7 in p, "cl_m7 is missing")
        self.failUnless(s7 in p, "s7 is missing")
        self.failUnless(cl_m8 in p, "cl_m8 is missing")
        self.failUnless(s8 in p, "s8 is missing")
        self.failUnless(cl_m9 in p, "cl_m9 is missing")
        self.failUnless(s9 in p, "s9 is missing")
        self.failUnless(cl_m10 in p, "cl_m10 is missing")
        self.failUnless(s10 in p, "s10 is missing")

        cp=s.copy()
        self.failUnless(isinstance(cp,SurfaceLoop), "copy returns is not an arc.")
        self.failUnless(not cp == s, "copy equals source")
        self.failUnless(cp.isColocated(s),"copy is not colocated with its source.")
        cc=cp.getSurfaces()
        self.failUnless(len(cc) == 10, "too many primitives in copy.")
        self.failUnless(not s1 in cc,"copy uses s1.")
        self.failUnless(not s2 in cc,"copy uses s2.")
        self.failUnless(not s3 in cc,"copy uses s3.")
        self.failUnless(not s4 in cc,"copy uses s4.")
        self.failUnless(not s5 in cc,"copy uses s5.")
        self.failUnless(not s6 in cc,"copy uses s6.")
        self.failUnless(not s7 in cc,"copy uses s7.")
        self.failUnless(not s8 in cc,"copy uses s8.")
        self.failUnless(not s9 in cc,"copy uses s9.")
        self.failUnless(not s10 in cc,"copy uses s10.")

        #=================================================================================
        p0_m=Point( 0, 0, 0,0.1)
        p1_m=Point(-10, 0, 0,0.1)
        p2_m=Point( 0,-10, 0,0.1)
        p3_m=Point(-10,-10, 0,0.1)
        p4_m=Point( 0, 0,-10,0.1)
        p5_m=Point(-10, 0,-10,0.1)
        p6_m=Point( 0,-10,-10,0.1)
        p7_m=Point(-10,-10,-10,0.1)
        q0_m=Point( -4, 0, -4,0.1)
        q1_m=Point( -6, 0, -4,0.1)
        q2_m=Point( -4,-10, -4,0.1)
        q3_m=Point( -6,-10, -4,0.1)
        q4_m=Point( -4, 0, -6,0.1)
        q5_m=Point( -6, 0, -6,0.1)
        q6_m=Point( -4,-10, -6,0.1)
        q7_m=Point( -6,-10, -6,0.1)
        l01_m=Line(p0_m,p1_m)
        l15_m=Line(p1_m,p5_m)
        l54_m=Line(p5_m,p4_m)
        l04_m=Line(p4_m,p0_m)
        l13_m=Line(p1_m,p3_m)
        l37_m=Line(p3_m,p7_m)
        l75_m=Line(p7_m,p5_m)
        l67_m=Line(p6_m,p7_m)
        l26_m=Line(p2_m,p6_m)
        l32_m=Line(p3_m,p2_m)
        l20_m=Line(p2_m,p0_m)
        l46_m=Line(p4_m,p6_m)
        m01_m=Line(q0_m,q1_m)
        m15_m=Line(q1_m,q5_m)
        m54_m=Line(q5_m,q4_m)
        m40_m=Line(q4_m,q0_m)
        m23_m=Line(q2_m,q3_m)
        m37_m=Line(q3_m,q7_m)
        m76_m=Line(q7_m,q6_m)
        m62_m=Line(q6_m,q2_m)
        m02_m=Line(q0_m,q2_m)
        m13_m=Line(q1_m,q3_m)
        m46_m=Line(q4_m,q6_m)
        m57_m=Line(q5_m,q7_m)
        cl_l1_m=CurveLoop(l01_m,l15_m,l54_m,l04_m)
        cl_m1_m=CurveLoop(m01_m,m15_m,m54_m,m40_m)
        s1_m=PlaneSurface(cl_l1_m,holes=[cl_m1_m])
        cl_l2_m=CurveLoop(-l15_m,l13_m,l37_m,l75_m)
        s2_m=PlaneSurface(cl_l2_m)
        cl_l3_m=CurveLoop(l32_m,-l37_m,l67_m,l26_m)
        cl_m3_m=CurveLoop(-m23_m,-m37_m,-m76_m,-m62_m)
        s3_m=PlaneSurface(cl_l3_m,holes=[cl_m3_m])
        cl_l4_m=CurveLoop(l20_m,-l26_m,l46_m,-l04_m)
        s4_m=PlaneSurface(cl_l4_m)
        cl_l5_m=CurveLoop(l32_m,l20_m,l01_m,l13_m)
        s5_m=PlaneSurface(-cl_l5_m)
        cl_l6_m=CurveLoop(l67_m,l75_m,l54_m,l46_m)
        s6_m=PlaneSurface(-cl_l6_m)
        cl_m7_m=CurveLoop(m13_m,m37_m,-m57_m,-m15_m)
        s7_m=PlaneSurface(cl_m7_m)
        cl_m8_m=CurveLoop(m57_m,m76_m,-m46_m,-m54_m)
        s8_m=PlaneSurface(cl_m8_m)
        cl_m9_m=CurveLoop(m46_m,m62_m,-m02_m,-m40_m)
        s9_m=PlaneSurface(cl_m9_m)
        cl_m10_m=CurveLoop(-m01_m,m02_m,m23_m,-m13_m)
        s10_m=PlaneSurface(cl_m10_m)
        s_m=SurfaceLoop(s1_m,s2_m,s3_m,s4_m,s5_m,s6_m,-s7_m,-s8_m,-s9_m,-s10_m)
        #_m=================================================================================
        dc=s.apply(Dilation(-1.))
        self.failUnless(dc.isColocated(s_m),"dilation is wrong.")
        cc=dc.getSurfaces()
        self.failUnless(len(cc) == 10, "too many surfaces in copy.")
        self.failUnless(not s1 in cc,"dilation uses s1.")
        self.failUnless(not s2 in cc,"dilation uses s2.")
        self.failUnless(not s3 in cc,"dilation uses s3.")
        self.failUnless(not s4 in cc,"dilation uses s4.")
        self.failUnless(not s5 in cc,"dilation uses s5.")
        self.failUnless(not s6 in cc,"dilation uses s6.")
        self.failUnless(not s7 in cc,"dilation uses s7.")
        self.failUnless(not s8 in cc,"dilation uses s8.")
        self.failUnless(not s9 in cc,"dilation uses s9.")
        self.failUnless(not s10 in cc,"dilation uses s10.")

        s.modifyBy(Dilation(-1.))
        self.failUnless(s.isColocated(s_m),"inplace dilation is wrong.")
        cc=s.getSurfaces()
        self.failUnless(len(cc) == 10, "too many primitives in modified object.")
        self.failUnless(s1 in cc,"dilation misses s1.")
        self.failUnless(cc[cc.index(s1)].hasSameOrientation(s1),"s1 in modified object has wrong orientation.")
        self.failUnless(s1.isColocated(s1_m),"s1 in modified object as wrong location.")
        self.failUnless(s2 in cc,"dilation misses s2.")
        self.failUnless(cc[cc.index(s2)].hasSameOrientation(s2),"s2 in modified object has wrong orientation.")
        self.failUnless(s2.isColocated(s2_m),"s2 in modified object as wrong location.")
        self.failUnless(s3 in cc,"dilation misses s3.")
        self.failUnless(cc[cc.index(s3)].hasSameOrientation(s3),"s3 in modified object has wrong orientation.")
        self.failUnless(s3.isColocated(s3_m),"s3 in modified object as wrong location.")
        self.failUnless(s4 in cc,"dilation misses s4.")
        self.failUnless(cc[cc.index(s4)].hasSameOrientation(s4),"s4 in modified object has wrong orientation.")
        self.failUnless(s4.isColocated(s4_m),"s4 in modified object as wrong location.")
        self.failUnless(s5 in cc,"dilation misses s5.")
        self.failUnless(cc[cc.index(s5)].hasSameOrientation(s5),"s5 in modified object has wrong orientation.")
        self.failUnless(s5.isColocated(s5_m),"s5 in modified object as wrong location.")
        self.failUnless(s6 in cc,"dilation misses s6.")
        self.failUnless(cc[cc.index(s6)].hasSameOrientation(s6),"s6 in modified object has wrong orientation.")
        self.failUnless(s6.isColocated(s6_m),"s6 in modified object as wrong location.")
        self.failUnless(s7 in cc,"dilation misses s7.")
        self.failUnless(cc[cc.index(s7)].hasSameOrientation(-s7),"s7 in modified object has wrong orientation.")
        self.failUnless(s7.isColocated(s7_m),"s7 in modified object as wrong location.")
        self.failUnless(s8 in cc,"dilation misses s8.")
        self.failUnless(cc[cc.index(s8)].hasSameOrientation(-s8),"s8 in modified object has wrong orientation.")
        self.failUnless(s8.isColocated(s8_m),"s8 in modified object as wrong location.")
        self.failUnless(s9 in cc,"dilation misses s9.")
        self.failUnless(cc[cc.index(s9)].hasSameOrientation(-s9),"s9 in modified object has wrong orientation.")
        self.failUnless(s9.isColocated(s9_m),"s9 in modified object as wrong location.")
        self.failUnless(s10 in cc,"dilation misses s10.")
        self.failUnless(cc[cc.index(s10)].hasSameOrientation(-s10),"s10 in modified object has wrong orientation.")
        self.failUnless(s10.isColocated(s10_m),"s10 in modified object as wrong location.")

   def test_ReverseSurfaceLoop(self):
        p0=Point( 0, 0, 0,0.1)
        p1=Point(10, 0, 0,0.1)
        p2=Point( 0,10, 0,0.1)
        p3=Point(10,10, 0,0.1)
        p4=Point( 0, 0,10,0.1)
        p5=Point(10, 0,10,0.1)
        p6=Point( 0,10,10,0.1)
        p7=Point(10,10,10,0.1)

        q0=Point( 4, 0, 4,0.1)
        q1=Point( 6, 0, 4,0.1)
        q2=Point( 4,10, 4,0.1)
        q3=Point( 6,10, 4,0.1)
        q4=Point( 4, 0, 6,0.1)
        q5=Point( 6, 0, 6,0.1)
        q6=Point( 4,10, 6,0.1)
        q7=Point( 6,10, 6,0.1)

        l01=Line(p0,p1)
        l15=Line(p1,p5)
        l54=Line(p5,p4)
        l04=Line(p4,p0)

        l13=Line(p1,p3)
        l37=Line(p3,p7)
        l75=Line(p7,p5)
        l67=Line(p6,p7)
        l26=Line(p2,p6)
        l32=Line(p3,p2)
        l20=Line(p2,p0)
        l46=Line(p4,p6)

        m01=Line(q0,q1)
        m15=Line(q1,q5)
        m54=Line(q5,q4)
        m40=Line(q4,q0)
        m23=Line(q2,q3)
        m37=Line(q3,q7)
        m76=Line(q7,q6)
        m62=Line(q6,q2)

        m02=Line(q0,q2)
        m13=Line(q1,q3)
        m46=Line(q4,q6)
        m57=Line(q5,q7)

        cl_l1=CurveLoop(l01,l15,l54,l04)
        cl_m1=CurveLoop(m01,m15,m54,m40)
        s1=PlaneSurface(cl_l1,holes=[cl_m1])
        s1_v=PlaneSurface(cl_l1)

        cl_l2=CurveLoop(-l15,l13,l37,l75)
        s2=PlaneSurface(cl_l2)

        cl_l3=CurveLoop(l32,-l37,l67,l26)
        cl_m3=CurveLoop(-m23,-m37,-m76,-m62)
        s3=PlaneSurface(cl_l3,holes=[cl_m3])
        s3_v=PlaneSurface(cl_l3)
     
     
        cl_l4=CurveLoop(l20,-l26,l46,-l04)
        s4=PlaneSurface(cl_l4)

        cl_l5=CurveLoop(l32,l20,l01,l13)
        s5=PlaneSurface(-cl_l5)

        cl_l6=CurveLoop(l67,l75,l54,l46)
        s6=PlaneSurface(-cl_l6)

        cl_m7=CurveLoop(m13,m37,-m57,-m15)
        s7=PlaneSurface(cl_m7)
        
        cl_m8=CurveLoop(m57,m76,-m46,-m54)
        s8=PlaneSurface(cl_m8)

        cl_m9=CurveLoop(m46,m62,-m02,-m40)
        s9=PlaneSurface(cl_m9)

        cl_m10=CurveLoop(-m01,m02,m23,-m13)
        s10=PlaneSurface(cl_m10)

        # self.failUnlessRaises(ValueError,SurfaceLoop,s1,s3)
        # self.failUnlessRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5)
        # self.failUnlessRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5,s5)

        
        CC0=SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)
        s=-CC0

        self.failUnless(s.hasSameOrientation(s),"has not same orientation like itself")
        self.failUnless(not s.hasSameOrientation(-s),"has same orientation like -itself")

        cc=s.getSurfaces()
        self.failUnless(len(cc) == 10, "too many curves.")
        self.failUnless(s1 in cc, "s1 is missing")
        self.failUnless(s2 in cc, "s2 is missing")
        self.failUnless(s3 in cc, "s3 is missing")
        self.failUnless(s4 in cc, "s4 is missing")
        self.failUnless(s5 in cc, "s5 is missing")
        self.failUnless(s6 in cc, "s6 is missing")
        self.failUnless(s7 in cc, "s7 is missing")
        self.failUnless(s8 in cc, "s8 is missing")
        self.failUnless(s9 in cc, "s9 is missing")
        self.failUnless(s10 in cc, "s10 is missing")

        self.failUnless(not s.isColocated(p4),"SurfaceLoop is colocated with point.")
        self.failUnless(s.isColocated(s),"SurfaceLoop is not colocated with its self.")
        self.failUnless(s.isColocated(-s),"SurfaceLoop is not colocated with its reverse.")
        self.failUnless(s.isColocated(SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)),"SurfaceLoop is not colocated with its copy.")
        self.failUnless(s.isColocated(SurfaceLoop(-s10, s1,s2,s3,s4,s5,s6,-s7,-s8,-s9)),"SurfaceLoop is not colocated with its copy with shifted points.")
        self.failUnless(s.isColocated(SurfaceLoop(s1, -s7, s2, -s9, s3, -s8, s4, -s10, s5,s6)),"SurfaceLoop is not colocated with its copy with shuffled points.")
        self.failUnless(not s.isColocated(SurfaceLoop(s1_v,s2,s3_v,s4,s5,s6)),"SurfaceLoop is colocated with different SurfaceLoop.")

        self.failUnless(len(s) == 10, "wrong length")

        s.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.failUnless(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.failUnless(p6.getLocalScale()==3., "p6 has wrong local scale.")
        self.failUnless(p7.getLocalScale()==3., "p7 has wrong local scale.")

        p=s.getPrimitives()
        self.failUnless(len(p) == 63, "too many primitives.")
        self.failUnless(s in p, "s is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")
        self.failUnless(p5 in p, "p5 is missing")
        self.failUnless(p6 in p, "p6 is missing")
        self.failUnless(p7 in p, "p7 is missing")
        self.failUnless(q0 in p, "q0 is missing")
        self.failUnless(q1 in p, "q1 is missing")
        self.failUnless(q2 in p, "q2 is missing")
        self.failUnless(q3 in p, "q3 is missing")
        self.failUnless(q4 in p, "q4 is missing")
        self.failUnless(q5 in p, "q5 is missing")
        self.failUnless(q6 in p, "q6 is missing")
        self.failUnless(q7 in p, "q7 is missing")
        self.failUnless(l01 in p, "l01 is missing")
        self.failUnless(l15 in p, "l15 is missing")
        self.failUnless(l54 in p, "l54 is missing")
        self.failUnless(l04 in p, "l04 is missing")
        self.failUnless(l13 in p, "l13 is missing")
        self.failUnless(l37 in p, "l37 is missing")
        self.failUnless(l75 in p, "l75 is missing")
        self.failUnless(l67 in p, "l67 is missing")
        self.failUnless(l26 in p, "l26 is missing")
        self.failUnless(l32 in p, "l32 is missing")
        self.failUnless(l20 in p, "l20 is missing")
        self.failUnless(l46 in p, "l46 is missing")
        self.failUnless(m01 in p, "m01 is missing")
        self.failUnless(m15 in p, "m15 is missing")
        self.failUnless(m54 in p, "m54 is missing")
        self.failUnless(m40 in p, "m40 is missing")
        self.failUnless(m23 in p, "m23 is missing")
        self.failUnless(m37 in p, "m37 is missing")
        self.failUnless(m76 in p, "m76 is missing")
        self.failUnless(m62 in p, "m62 is missing")
        self.failUnless(m02 in p, "m02 is missing")
        self.failUnless(m13 in p, "m13 is missing")
        self.failUnless(m46 in p, "m46 is missing")
        self.failUnless(m57 in p, "m57 is missing")
        self.failUnless(cl_l1 in p, "cl_l1 is missing")
        self.failUnless(cl_m1 in p, "cl_m1 is missing")
        self.failUnless(s1 in p, "s1 is missing")
        self.failUnless(cl_l2 in p, "cl_l2 is missing")
        self.failUnless(s2 in p, "s2 is missing")
        self.failUnless(cl_l3 in p, "cl_l3 is missing")
        self.failUnless(cl_m3 in p, "cl_m3 is missing")
        self.failUnless(s3 in p, "s3 is missing")
        self.failUnless(cl_l4 in p, "cl_l4 is missing")
        self.failUnless(s4 in p, "s4 is missing")
        self.failUnless(cl_l5 in p, "cl_l5 is missing")
        self.failUnless(s5 in p, "s5 is missing")
        self.failUnless(cl_l6 in p, "cl_l6 is missing")
        self.failUnless(s6 in p, "s6 is missing")
        self.failUnless(cl_m7 in p, "cl_m7 is missing")
        self.failUnless(s7 in p, "s7 is missing")
        self.failUnless(cl_m8 in p, "cl_m8 is missing")
        self.failUnless(s8 in p, "s8 is missing")
        self.failUnless(cl_m9 in p, "cl_m9 is missing")
        self.failUnless(s9 in p, "s9 is missing")
        self.failUnless(cl_m10 in p, "cl_m10 is missing")
        self.failUnless(s10 in p, "s10 is missing")

        cp=s.copy()
        self.failUnless(isinstance(cp,ReverseSurfaceLoop), "copy returns is not ReverseSurfaceLoop.")
        self.failUnless(not cp == s, "copy equals source")
        self.failUnless(cp.isColocated(s),"copy is not colocated with its source.")
        cc=cp.getSurfaces()
        self.failUnless(len(cc) == 10, "too many primitives in copy.")
        self.failUnless(not s1 in cc,"copy uses s1.")
        self.failUnless(not s2 in cc,"copy uses s2.")
        self.failUnless(not s3 in cc,"copy uses s3.")
        self.failUnless(not s4 in cc,"copy uses s4.")
        self.failUnless(not s5 in cc,"copy uses s5.")
        self.failUnless(not s6 in cc,"copy uses s6.")
        self.failUnless(not s7 in cc,"copy uses s7.")
        self.failUnless(not s8 in cc,"copy uses s8.")
        self.failUnless(not s9 in cc,"copy uses s9.")
        self.failUnless(not s10 in cc,"copy uses s10.")

        #=================================================================================
        p0_m=Point( 0, 0, 0,0.1)
        p1_m=Point(-10, 0, 0,0.1)
        p2_m=Point( 0,-10, 0,0.1)
        p3_m=Point(-10,-10, 0,0.1)
        p4_m=Point( 0, 0,-10,0.1)
        p5_m=Point(-10, 0,-10,0.1)
        p6_m=Point( 0,-10,-10,0.1)
        p7_m=Point(-10,-10,-10,0.1)
        q0_m=Point( -4, 0, -4,0.1)
        q1_m=Point( -6, 0, -4,0.1)
        q2_m=Point( -4,-10, -4,0.1)
        q3_m=Point( -6,-10, -4,0.1)
        q4_m=Point( -4, 0, -6,0.1)
        q5_m=Point( -6, 0, -6,0.1)
        q6_m=Point( -4,-10, -6,0.1)
        q7_m=Point( -6,-10, -6,0.1)
        l01_m=Line(p0_m,p1_m)
        l15_m=Line(p1_m,p5_m)
        l54_m=Line(p5_m,p4_m)
        l04_m=Line(p4_m,p0_m)
        l13_m=Line(p1_m,p3_m)
        l37_m=Line(p3_m,p7_m)
        l75_m=Line(p7_m,p5_m)
        l67_m=Line(p6_m,p7_m)
        l26_m=Line(p2_m,p6_m)
        l32_m=Line(p3_m,p2_m)
        l20_m=Line(p2_m,p0_m)
        l46_m=Line(p4_m,p6_m)
        m01_m=Line(q0_m,q1_m)
        m15_m=Line(q1_m,q5_m)
        m54_m=Line(q5_m,q4_m)
        m40_m=Line(q4_m,q0_m)
        m23_m=Line(q2_m,q3_m)
        m37_m=Line(q3_m,q7_m)
        m76_m=Line(q7_m,q6_m)
        m62_m=Line(q6_m,q2_m)
        m02_m=Line(q0_m,q2_m)
        m13_m=Line(q1_m,q3_m)
        m46_m=Line(q4_m,q6_m)
        m57_m=Line(q5_m,q7_m)
        cl_l1_m=CurveLoop(l01_m,l15_m,l54_m,l04_m)
        cl_m1_m=CurveLoop(m01_m,m15_m,m54_m,m40_m)
        s1_m=PlaneSurface(cl_l1_m,holes=[cl_m1_m])
        cl_l2_m=CurveLoop(-l15_m,l13_m,l37_m,l75_m)
        s2_m=PlaneSurface(cl_l2_m)
        cl_l3_m=CurveLoop(l32_m,-l37_m,l67_m,l26_m)
        cl_m3_m=CurveLoop(-m23_m,-m37_m,-m76_m,-m62_m)
        s3_m=PlaneSurface(cl_l3_m,holes=[cl_m3_m])
        cl_l4_m=CurveLoop(l20_m,-l26_m,l46_m,-l04_m)
        s4_m=PlaneSurface(cl_l4_m)
        cl_l5_m=CurveLoop(l32_m,l20_m,l01_m,l13_m)
        s5_m=PlaneSurface(-cl_l5_m)
        cl_l6_m=CurveLoop(l67_m,l75_m,l54_m,l46_m)
        s6_m=PlaneSurface(-cl_l6_m)
        cl_m7_m=CurveLoop(m13_m,m37_m,-m57_m,-m15_m)
        s7_m=PlaneSurface(cl_m7_m)
        cl_m8_m=CurveLoop(m57_m,m76_m,-m46_m,-m54_m)
        s8_m=PlaneSurface(cl_m8_m)
        cl_m9_m=CurveLoop(m46_m,m62_m,-m02_m,-m40_m)
        s9_m=PlaneSurface(cl_m9_m)
        cl_m10_m=CurveLoop(-m01_m,m02_m,m23_m,-m13_m)
        s10_m=PlaneSurface(cl_m10_m)
        s_m=SurfaceLoop(s1_m,s2_m,s3_m,s4_m,s5_m,s6_m,-s7_m,-s8_m,-s9_m,-s10_m)
        #_m=================================================================================
        dc=s.apply(Dilation(-1.))
        self.failUnless(dc.isColocated(s_m),"dilation is wrong.")
        cc=dc.getSurfaces()
        self.failUnless(len(cc) == 10, "too many surfaces in copy.")
        self.failUnless(not s1 in cc,"dilation uses s1.")
        self.failUnless(not s2 in cc,"dilation uses s2.")
        self.failUnless(not s3 in cc,"dilation uses s3.")
        self.failUnless(not s4 in cc,"dilation uses s4.")
        self.failUnless(not s5 in cc,"dilation uses s5.")
        self.failUnless(not s6 in cc,"dilation uses s6.")
        self.failUnless(not s7 in cc,"dilation uses s7.")
        self.failUnless(not s8 in cc,"dilation uses s8.")
        self.failUnless(not s9 in cc,"dilation uses s9.")
        self.failUnless(not s10 in cc,"dilation uses s10.")

        s.modifyBy(Dilation(-1.))
        self.failUnless(s.isColocated(s_m),"inplace dilation is wrong.")
        cc=s.getSurfaces()
        self.failUnless(len(cc) == 10, "too many primitives in modified object.")
        self.failUnless(s1 in cc,"dilation misses s1.")
        self.failUnless(cc[cc.index(s1)].hasSameOrientation(-s1),"s1 in modified object has wrong orientation.")
        self.failUnless(s1.isColocated(s1_m),"s1 in modified object as wrong location.")
        self.failUnless(s2 in cc,"dilation misses s2.")
        self.failUnless(cc[cc.index(s2)].hasSameOrientation(-s2),"s2 in modified object has wrong orientation.")
        self.failUnless(s2.isColocated(s2_m),"s2 in modified object as wrong location.")
        self.failUnless(s3 in cc,"dilation misses s3.")
        self.failUnless(cc[cc.index(s3)].hasSameOrientation(-s3),"s3 in modified object has wrong orientation.")
        self.failUnless(s3.isColocated(s3_m),"s3 in modified object as wrong location.")
        self.failUnless(s4 in cc,"dilation misses s4.")
        self.failUnless(cc[cc.index(s4)].hasSameOrientation(-s4),"s4 in modified object has wrong orientation.")
        self.failUnless(s4.isColocated(s4_m),"s4 in modified object as wrong location.")
        self.failUnless(s5 in cc,"dilation misses s5.")
        self.failUnless(cc[cc.index(s5)].hasSameOrientation(-s5),"s5 in modified object has wrong orientation.")
        self.failUnless(s5.isColocated(s5_m),"s5 in modified object as wrong location.")
        self.failUnless(s6 in cc,"dilation misses s6.")
        self.failUnless(cc[cc.index(s6)].hasSameOrientation(-s6),"s6 in modified object has wrong orientation.")
        self.failUnless(s6.isColocated(s6_m),"s6 in modified object as wrong location.")
        self.failUnless(s7 in cc,"dilation misses s7.")
        self.failUnless(cc[cc.index(s7)].hasSameOrientation(s7),"s7 in modified object has wrong orientation.")
        self.failUnless(s7.isColocated(s7_m),"s7 in modified object as wrong location.")
        self.failUnless(s8 in cc,"dilation misses s8.")
        self.failUnless(cc[cc.index(s8)].hasSameOrientation(s8),"s8 in modified object has wrong orientation.")
        self.failUnless(s8.isColocated(s8_m),"s8 in modified object as wrong location.")
        self.failUnless(s9 in cc,"dilation misses s9.")
        self.failUnless(cc[cc.index(s9)].hasSameOrientation(s9),"s9 in modified object has wrong orientation.")
        self.failUnless(s9.isColocated(s9_m),"s9 in modified object as wrong location.")
        self.failUnless(s10 in cc,"dilation misses s10.")
        self.failUnless(cc[cc.index(s10)].hasSameOrientation(s10),"s10 in modified object has wrong orientation.")
        self.failUnless(s10.isColocated(s10_m),"s10 in modified object as wrong location.")

   def test_Volume(self):
        p0=Point(-2,-2,-2,0.1)
        p1=Point(2,-2,-2,0.1)
        p2=Point(-2,2,-2,0.1)
        p3=Point(2,2,-2,0.1)
        p4=Point(-2,-2,2,0.1)
        p5=Point(2,-2,2,0.1)
        p6=Point(-2,2,2,0.1)
        p7=Point(2,2,2,0.1)
        l01=Line(p0,p1)
        l15=Line(p1,p5)
        l54=Line(p5,p4)
        l40=Line(p4,p0)
        l23=Line(p2,p3)
        l37=Line(p3,p7)
        l76=Line(p7,p6)
        l62=Line(p6,p2)
        l13=Line(p1,p3)
        l57=Line(p5,p7)
        l02=Line(p0,p2)
        l46=Line(p4,p6)
        cl1=CurveLoop(l01,l15,l54,l40)
        s1=PlaneSurface(cl1)
        cl2=CurveLoop(l23,l37,l76,l62)
        s2=PlaneSurface(-cl2)
        cl3=CurveLoop(l13,l37,-l57,-l15)
        s3=PlaneSurface(cl3)
        cl4=CurveLoop(l46,l62,-l02,-l40)
        s4=PlaneSurface(-cl4)
        cl5=CurveLoop(-l01,l02,l23,-l13)
        s5=PlaneSurface(-cl5)
        cl6=CurveLoop(-l54,l57,l76,-l46)
        s6=PlaneSurface(-cl6)
        s_out=SurfaceLoop(s1,s2,s3,s4,s5,s6)
 
        p0_i=Point(-1,-1,-1,0.1)
        p1_i=Point(1,-1,-1,0.1)
        p2_i=Point(-1,1,-1,0.1)
        p3_i=Point(1,1,-1,0.1)
        p4_i=Point(-1,-1,1,0.1)
        p5_i=Point(1,-1,1,0.1)
        p6_i=Point(-1,1,1,0.1)
        p7_i=Point(1,1,1,0.1)
        l01_i=Line(p0_i,p1_i)
        l15_i=Line(p1_i,p5_i)
        l54_i=Line(p5_i,p4_i)
        l40_i=Line(p4_i,p0_i)
        l23_i=Line(p2_i,p3_i)
        l37_i=Line(p3_i,p7_i)
        l76_i=Line(p7_i,p6_i)
        l62_i=Line(p6_i,p2_i)
        l13_i=Line(p1_i,p3_i)
        l57_i=Line(p5_i,p7_i)
        l02_i=Line(p0_i,p2_i)
        l46_i=Line(p4_i,p6_i)
        cl1_i=CurveLoop(l01_i,l15_i,l54_i,l40_i)
        s1_i=PlaneSurface(cl1_i)
        cl2_i=CurveLoop(l23_i,l37_i,l76_i,l62_i)
        s2_i=PlaneSurface(-cl2_i)
        cl3_i=CurveLoop(l13_i,l37_i,-l57_i,-l15_i)
        s3_i=PlaneSurface(cl3_i)
        cl4_i=CurveLoop(l46_i,l62_i,-l02_i,-l40_i)
        s4_i=PlaneSurface(-cl4_i)
        cl5_i=CurveLoop(-l01_i,l02_i,l23_i,-l13_i)
        s5_i=PlaneSurface(-cl5_i)
        cl6_i=CurveLoop(-l54_i,l57_i,l76_i,-l46_i)
        s6_i=PlaneSurface(-cl6_i)
        s_inner=SurfaceLoop(s1_i,s2_i,s3_i,s4_i,s5_i,s6_i)


        self.failUnlessRaises(TypeError,Volume,s1)
        self.failUnlessRaises(TypeError,Volume,s_out,[s1])
        v=Volume(s_out,holes=[s_inner])

        self.failUnlessRaises(NotImplementedError,v.__neg__)
        
        cl2=v.getSurfaceLoop()
        self.failUnless(s_out == cl2, " wrong boundary loops")
        self.failUnless(s_out.hasSameOrientation(cl2),"cl has incorrect orientation.")

        hs=v.getHoles()
        self.failUnless(len(hs) == 1, "one holes expected.")
        self.failUnless(s_inner==hs[0], "h is not defined as hole.")
        self.failUnless(hs[0].hasSameOrientation(s_inner),"hs has incorrect orientation.")

        cc=v.getBoundary()
        self.failUnless(len(cc) == 12, "too many curves.")
        self.failUnless(s1 in cc, "s1 is missing")
        self.failUnless(s2 in cc, "s2 is missing")
        self.failUnless(s3 in cc, "s3 is missing")
        self.failUnless(s4 in cc, "s4 is missing")
        self.failUnless(s5 in cc, "s5 is missing")
        self.failUnless(s6 in cc, "s6 is missing")
        self.failUnless(s1_i in cc, "s1_i is missing")
        self.failUnless(s2_i in cc, "s2_i is missing")
        self.failUnless(s3_i in cc, "s3_i is missing")
        self.failUnless(s4_i in cc, "s4_i is missing")
        self.failUnless(s5_i in cc, "s5_i is missing")
        self.failUnless(s6_i in cc, "s6_i is missing")

        self.failUnless(not v.isColocated(p4),"Volume is colocated with point.")
        self.failUnless(v.isColocated(v),"Volume is not colocated with its self.")
        self.failUnless(v.isColocated(Volume(s_out,holes=[s_inner])),"Volume is not colocated with its copy.")
        self.failUnless(not v.isColocated(Volume(s_out,holes=[0.3*s_inner])),"Volume is not colocated with volume of modifies volume .")

        v.setLocalScale(3.)
        self.failUnless(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.failUnless(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.failUnless(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.failUnless(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.failUnless(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.failUnless(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.failUnless(p6.getLocalScale()==3., "p6 has wrong local scale.")
        self.failUnless(p7.getLocalScale()==3., "p7 has wrong local scale.")
        self.failUnless(p0_i.getLocalScale()==3., "p0_i has wrong local scale.")
        self.failUnless(p1_i.getLocalScale()==3., "p1_i has wrong local scale.")
        self.failUnless(p2_i.getLocalScale()==3., "p2_i has wrong local scale.")
        self.failUnless(p3_i.getLocalScale()==3., "p3_i has wrong local scale.")
        self.failUnless(p4_i.getLocalScale()==3., "p4_i has wrong local scale.")
        self.failUnless(p5_i.getLocalScale()==3., "p5_i has wrong local scale.")
        self.failUnless(p6_i.getLocalScale()==3., "p6_i has wrong local scale.")
        self.failUnless(p7_i.getLocalScale()==3., "p7_i has wrong local scale.")

        p=v.getPrimitives()
        self.failUnless(len(p) == 67, "too many primitives.")
        self.failUnless(v in p, "v is missing")
        self.failUnless(p0 in p, "p0 is missing")
        self.failUnless(p1 in p, "p1 is missing")
        self.failUnless(p2 in p, "p2 is missing")
        self.failUnless(p3 in p, "p3 is missing")
        self.failUnless(p4 in p, "p4 is missing")
        self.failUnless(p5 in p, "p5 is missing")
        self.failUnless(p6 in p, "p6 is missing")
        self.failUnless(p7 in p, "p7 is missing")
        self.failUnless(l01 in p, "l01 is missing")
        self.failUnless(l15 in p, "l15 is missing")
        self.failUnless(l54 in p, "l54 is missing")
        self.failUnless(l40 in p, "l40 is missing")
        self.failUnless(l23 in p, "l23 is missing")
        self.failUnless(l37 in p, "l37 is missing")
        self.failUnless(l76 in p, "l76 is missing")
        self.failUnless(l62 in p, "l62 is missing")
        self.failUnless(l13 in p, "l13 is missing")
        self.failUnless(l57 in p, "l57 is missing")
        self.failUnless(l02 in p, "l02 is missing")
        self.failUnless(l46 in p, "l46 is missing")
        self.failUnless(cl1 in p, "cl1 is missing")
        self.failUnless(s1 in p, "s1 is missing")
        self.failUnless(cl2 in p, "cl2 is missing")
        self.failUnless(s2 in p, "s2 is missing")
        self.failUnless(cl3 in p, "cl3  is missing")
        self.failUnless(s3 in p, "s3  is missing")
        self.failUnless(cl4 in p, "cl4  is missing")
        self.failUnless(s4 in p, "s4  is missing")
        self.failUnless(cl5 in p, "cl5  is missing")
        self.failUnless(s5 in p, "s5  is missing")
        self.failUnless(cl6 in p, "cl6  is missing")
        self.failUnless(s6 in p, "s6  is missing")
        self.failUnless(s_out in p, "s_out is missing")
        self.failUnless(p0_i in p, "p0_i is missing")
        self.failUnless(p1_i in p, "p1_i is missing")
        self.failUnless(p2_i in p, "p2_i is missing")
        self.failUnless(p3_i in p, "p3_i is missing")
        self.failUnless(p4_i in p, "p4_i is missing")
        self.failUnless(p5_i in p, "p5_i is missing")
        self.failUnless(p6_i in p, "p6_i is missing")
        self.failUnless(p7_i in p, "p7_i is missing")
        self.failUnless(l01_i in p, "l01_i is missing")
        self.failUnless(l15_i in p, "l15_i is missing")
        self.failUnless(l54_i in p, "l54_i is missing")
        self.failUnless(l40_i in p, "l40_i is missing")
        self.failUnless(l23_i in p, "l23_i is missing")
        self.failUnless(l37_i in p, "l37_i is missing")
        self.failUnless(l76_i in p, "l76_i is missing")
        self.failUnless(l62_i in p, "l62_i is missing")
        self.failUnless(l13_i in p, "l13_i is missing")
        self.failUnless(l57_i in p, "l57_i is missing")
        self.failUnless(l02_i in p, "l02_i is missing")
        self.failUnless(l46_i in p, "l46_i is missing")
        self.failUnless(cl1_i in p, "cl1_i is missing")
        self.failUnless(s1_i in p, "s1_i is missing")
        self.failUnless(cl2_i in p, "cl2_i is missing")
        self.failUnless(s2_i in p, "s2_i is missing")
        self.failUnless(cl3_i in p, "cl3_i  is missing")
        self.failUnless(s3_i in p, "s3_i  is missing")
        self.failUnless(cl4_i in p, "cl4_i  is missing")
        self.failUnless(s4_i in p, "s4_i  is missing")
        self.failUnless(cl5_i in p, "cl5_i  is missing")
        self.failUnless(s5_i in p, "s5_i  is missing")
        self.failUnless(cl6_i in p, "cl6_i  is missing")
        self.failUnless(s6_i in p, "s6_i  is missing")
        self.failUnless(s_inner in p, "s_inner  is missing")

        cp=v.copy()
        self.failUnless(isinstance(cp,Volume), "copy returns is not an arc.")
        self.failUnless(not cp == v, "copy equals source")
        self.failUnless(cp.isColocated(v),"copy is not colocated with its source.")
        cl2=cp.getSurfaceLoop()
        self.failUnless(not s_out == cl2, "copy uses s_out")
        hs=cp.getHoles()
        self.failUnless(len(hs) == 1, "copy: one holes expected.")
        self.failUnless(not s_inner==hs[0], "copy: uses s_inner.")
        cc=cp.getBoundary()
        self.failUnless(len(cc) == 12, "too many primitives in copy.")
        self.failUnless(not s1 in cc,"copy uses s1.")
        self.failUnless(not s2 in cc,"copy uses s2.")
        self.failUnless(not s3 in cc,"copy uses s3.")
        self.failUnless(not s4 in cc,"copy uses s4.")
        self.failUnless(not s5 in cc,"copy uses s5.")
        self.failUnless(not s6 in cc,"copy uses s6.")
        self.failUnless(not s1_i in cc,"copy uses s1_i.")
        self.failUnless(not s2_i in cc,"copy uses s2_i.")
        self.failUnless(not s3_i in cc,"copy uses s3_i.")
        self.failUnless(not s4_i in cc,"copy uses s4_i.")
        self.failUnless(not s5_i in cc,"copy uses s5_i.")
        self.failUnless(not s6_i in cc,"copy uses s6_i.")

        #=================================================================================
        p0_m=Point(2,2,2,0.1)
        p1_m=Point(-2,2,2,0.1)
        p2_m=Point(2,-2,2,0.1)
        p3_m=Point(-2,-2,2,0.1)
        p4_m=Point(2,2,-2,0.1)
        p5_m=Point(-2,2,-2,0.1)
        p6_m=Point(2,-2,-2,0.1)
        p7_m=Point(-2,-2,-2,0.1)
        l01_m=Line(p0_m,p1_m)
        l15_m=Line(p1_m,p5_m)
        l54_m=Line(p5_m,p4_m)
        l40_m=Line(p4_m,p0_m)
        l23_m=Line(p2_m,p3_m)
        l37_m=Line(p3_m,p7_m)
        l76_m=Line(p7_m,p6_m)
        l62_m=Line(p6_m,p2_m)
        l13_m=Line(p1_m,p3_m)
        l57_m=Line(p5_m,p7_m)
        l02_m=Line(p0_m,p2_m)
        l46_m=Line(p4_m,p6_m)
        cl1_m=CurveLoop(l01_m,l15_m,l54_m,l40_m)
        s1_m=PlaneSurface(cl1_m)
        cl2_m=CurveLoop(l23_m,l37_m,l76_m,l62_m)
        s2_m=PlaneSurface(-cl2_m)
        cl3_m=CurveLoop(l13_m,l37_m,-l57_m,-l15_m)
        s3_m=PlaneSurface(cl3_m)
        cl4_m=CurveLoop(l46_m,l62_m,-l02_m,-l40_m)
        s4_m=PlaneSurface(-cl4_m)
        cl5_m=CurveLoop(-l01_m,l02_m,l23_m,-l13_m)
        s5_m=PlaneSurface(-cl5_m)
        cl6_m=CurveLoop(-l54_m,l57_m,l76_m,-l46_m)
        s6_m=PlaneSurface(-cl6_m)
        s_out_m=SurfaceLoop(s1_m,s2_m,s3_m,s4_m,s5_m,s6_m)
        p0_i_m=Point(1,1,1,0.1)
        p1_i_m=Point(-1,1,1,0.1)
        p2_i_m=Point(1,-1,1,0.1)
        p3_i_m=Point(-1,-1,1,0.1)
        p4_i_m=Point(1,1,-1,0.1)
        p5_i_m=Point(-1,1,-1,0.1)
        p6_i_m=Point(1,-1,-1,0.1)
        p7_i_m=Point(-1,-1,-1,0.1)
        l01_i_m=Line(p0_i_m,p1_i_m)
        l15_i_m=Line(p1_i_m,p5_i_m)
        l54_i_m=Line(p5_i_m,p4_i_m)
        l40_i_m=Line(p4_i_m,p0_i_m)
        l23_i_m=Line(p2_i_m,p3_i_m)
        l37_i_m=Line(p3_i_m,p7_i_m)
        l76_i_m=Line(p7_i_m,p6_i_m)
        l62_i_m=Line(p6_i_m,p2_i_m)
        l13_i_m=Line(p1_i_m,p3_i_m)
        l57_i_m=Line(p5_i_m,p7_i_m)
        l02_i_m=Line(p0_i_m,p2_i_m)
        l46_i_m=Line(p4_i_m,p6_i_m)
        cl1_i_m=CurveLoop(l01_i_m,l15_i_m,l54_i_m,l40_i_m)
        s1_i_m=PlaneSurface(cl1_i_m)
        cl2_i_m=CurveLoop(l23_i_m,l37_i_m,l76_i_m,l62_i_m)
        s2_i_m=PlaneSurface(-cl2_i_m)
        cl3_i_m=CurveLoop(l13_i_m,l37_i_m,-l57_i_m,-l15_i_m)
        s3_i_m=PlaneSurface(cl3_i_m)
        cl4_i_m=CurveLoop(l46_i_m,l62_i_m,-l02_i_m,-l40_i_m)
        s4_i_m=PlaneSurface(-cl4_i_m)
        cl5_i_m=CurveLoop(-l01_i_m,l02_i_m,l23_i_m,-l13_i_m)
        s5_i_m=PlaneSurface(-cl5_i_m)
        cl6_i_m=CurveLoop(-l54_i_m,l57_i_m,l76_i_m,-l46_i_m)
        s6_i_m=PlaneSurface(-cl6_i_m)
        s_inner_m=SurfaceLoop(s1_i_m,s2_i_m,s3_i_m,s4_i_m,s5_i_m,s6_i_m)
        v_m=Volume(s_out_m,holes=[s_inner_m])
        #===========================================================================================================
        dc=v.apply(Dilation(-1.))
        self.failUnless(dc.isColocated(v),"dilation is wrong.")
        cl2=cp.getSurfaceLoop()
        self.failUnless(not s_out == cl2, "copy uses s_out")
        hs=cp.getHoles()
        self.failUnless(len(hs) == 1, "copy: one holes expected.")
        self.failUnless(not s_inner==hs[0], "copy: uses s_inner.")
        cc=cp.getBoundary()
        self.failUnless(len(cc) == 12, "too many primitives in copy.")
        self.failUnless(not s1 in cc,"copy uses s1.")
        self.failUnless(not s2 in cc,"copy uses s2.")
        self.failUnless(not s3 in cc,"copy uses s3.")
        self.failUnless(not s4 in cc,"copy uses s4.")
        self.failUnless(not s5 in cc,"copy uses s5.")
        self.failUnless(not s6 in cc,"copy uses s6.")
        self.failUnless(not s1_i in cc,"copy uses s1_i.")
        self.failUnless(not s2_i in cc,"copy uses s2_i.")
        self.failUnless(not s3_i in cc,"copy uses s3_i.")
        self.failUnless(not s4_i in cc,"copy uses s4_i.")
        self.failUnless(not s5_i in cc,"copy uses s5_i.")
        self.failUnless(not s6_i in cc,"copy uses s6_i.")

        v.modifyBy(Dilation(-1.))
        self.failUnless(v.isColocated(v_m),"inplace dilation is wrong.")
        cl2=v.getSurfaceLoop()
        self.failUnless(s_out == cl2, "inplace dilation does not s_out")
        hs=v.getHoles()
        self.failUnless(len(hs) == 1, "inplace dilation: one holes expected.")
        self.failUnless(s_inner==hs[0], "inplace dilation does not s_inner as hole")
        cc=v.getBoundary()
        self.failUnless(len(cc) == 12, "too many primitives in copy.")
        self.failUnless(s1 in cc,"inplace dilation does not use s1.")
        self.failUnless(cc[cc.index(s1)].hasSameOrientation(s1),"s1 in modified object has wrong orientation.")
        self.failUnless(s1.isColocated(s1_m),"s1 in dilated object as wrong location.")
        self.failUnless(s2 in cc,"inplace dilation does not use s2.")
        self.failUnless(cc[cc.index(s2)].hasSameOrientation(s2),"s2 in modified object has wrong orientation.")
        self.failUnless(s2.isColocated(s2_m),"s2 in dilated object as wrong location.")
        self.failUnless(s3 in cc,"inplace dilation does not use s3.")
        self.failUnless(cc[cc.index(s3)].hasSameOrientation(s3),"s3 in modified object has wrong orientation.")
        self.failUnless(s3.isColocated(s3_m),"s3 in dilated object as wrong location.")
        self.failUnless(s4 in cc,"inplace dilation does not use s4.")
        self.failUnless(cc[cc.index(s4)].hasSameOrientation(s4),"s4 in modified object has wrong orientation.")
        self.failUnless(s4.isColocated(s4_m),"s4 in dilated object as wrong location.")
        self.failUnless(s5 in cc,"inplace dilation does not use s5.")
        self.failUnless(cc[cc.index(s5)].hasSameOrientation(s5),"s5 in modified object has wrong orientation.")
        self.failUnless(s5.isColocated(s5_m),"s5 in dilated object as wrong location.")
        self.failUnless(s6 in cc,"inplace dilation does not use s6.")
        self.failUnless(cc[cc.index(s6)].hasSameOrientation(s6),"s6 in modified object has wrong orientation.")
        self.failUnless(s6.isColocated(s6_m),"s6 in dilated object as wrong location.")
        self.failUnless(s1_i in cc,"inplace dilation does not use s1_i.")
        self.failUnless(cc[cc.index(s1_i)].hasSameOrientation(s1_i),"s1_i in modified object has wrong orientation.")
        self.failUnless(s1_i.isColocated(s1_i_m),"s1_i in dilated object as wrong location.")
        self.failUnless(s2_i in cc,"inplace dilation does not use s2_i.")
        self.failUnless(cc[cc.index(s2_i)].hasSameOrientation(s2_i),"s2_i in modified object has wrong orientation.")
        self.failUnless(s2_i.isColocated(s2_i_m),"s2_i in dilated object as wrong location.")
        self.failUnless(s3_i in cc,"inplace dilation does not use s3_i.")
        self.failUnless(cc[cc.index(s3_i)].hasSameOrientation(s3_i),"s3_i in modified object has wrong orientation.")
        self.failUnless(s3_i.isColocated(s3_i_m),"s3_i in dilated object as wrong location.")
        self.failUnless(s4_i in cc,"inplace dilation does not use s4_i.")
        self.failUnless(cc[cc.index(s4_i)].hasSameOrientation(s4_i),"s4_i in modified object has wrong orientation.")
        self.failUnless(s4_i.isColocated(s4_i_m),"s4_i in dilated object as wrong location.")
        self.failUnless(s5_i in cc,"inplace dilation does not use s5_i.")
        self.failUnless(cc[cc.index(s5_i)].hasSameOrientation(s5_i),"s5_i in modified object has wrong orientation.")
        self.failUnless(s5_i.isColocated(s5_i_m),"s5_i in dilated object as wrong location.")
        self.failUnless(s6_i in cc,"inplace dilation does not use s6_i.")
        self.failUnless(cc[cc.index(s6_i)].hasSameOrientation(s6_i),"s6_i in modified object has wrong orientation.")
        self.failUnless(s6_i.isColocated(s6_i_m),"s6_i in dilated object as wrong location.")

   def test_PropertySet0D(self):
       p0=Point(1.,2.,3.,local_scale=9.)
       p1=Point(0.,0.,0.,local_scale=9.)
       p3=Point(8.,6,6,local_scale=9.)
       p4=Point(8.,6,-6,local_scale=9.)
       
       # create property set with dim:
       self.failUnlessRaises(TypeError,PropertySet,"test0",p0, Line(p0,p1))

       ps=PropertySet("test0",p0, p1)
       self.failUnless(ps.getManifoldClass() == Point, "wrong manifold")
       self.failUnless(ps.getDim() == 0, "wrong dimension")
       self.failUnless(ps.getName() == "test0", "wrong name")
       ps.setName("test1")
       self.failUnless(ps.getName() == "test1", "wrong new name")

       self.failUnless(ps.getTag() == 7, "wrong tag")
       
       self.failUnlessRaises(TypeError,ps.addItem, Line(p0,p1))
       ps.addItem(p3,p4)
       pp=ps.getItems()
       self.failUnless(len(pp) == 4, "wrong number of items")
       self.failUnless(p0 in pp, "p0 missing in items.")
       self.failUnless(p1 in pp, "p1 missing in items.")
       self.failUnless(p3 in pp, "p3 missing in items.")
       self.failUnless(p4 in pp, "p4 missing in items.")

       pp=ps.getPrimitives()
       self.failUnless(len(pp) == 5, "wrong number of items")
       self.failUnless(ps in pp, "ps missing in items.")
       self.failUnless(p0 in pp, "p0 missing in items.")
       self.failUnless(p1 in pp, "p1 missing in items.")
       self.failUnless(p3 in pp, "p3 missing in items.")
       self.failUnless(p4 in pp, "p4 missing in items.")
       
       ps.clearItems()
       self.failUnless(len(ps.getItems()) == 0, "cleaning items failed.")
       
       
   def test_PropertySet1D(self):
       p0=Point(1.,2.,3.,local_scale=9.)
       p1=Point(0.,0.,0.,local_scale=9.)
       p3=Point(8.,6,6,local_scale=9.)
       p4=Point(8.,6,-6,local_scale=9.)
       
       l0=Line(p0,p1)
       l1=Arc(p3,p1,p4)
       l2=Line(p4,p0)
       # create property set with dim:
       self.failUnlessRaises(TypeError,PropertySet,"test0",l0, p0)

       ps=PropertySet("test0", l0, l1)
       self.failUnless(ps.getManifoldClass() == Manifold1D, "wrong manifold")
       self.failUnless(ps.getDim() == 1, "wrong dimension")
       self.failUnless(ps.getName() == "test0", "wrong name")
       ps.setName("test1")
       self.failUnless(ps.getName() == "test1", "wrong new name")

       self.failUnless(ps.getTag() == 9, "wrong tag")
       
       self.failUnlessRaises(TypeError,ps.addItem, p0)
       ps.addItem(l2)
       pp=ps.getItems()
       self.failUnless(len(pp) == 3, "wrong number of items")
       self.failUnless(l0 in pp, "l0 missing in items.")
       self.failUnless(l1 in pp, "l1 missing in items.")
       self.failUnless(l2 in pp, "l2 missing in items.")

       pp=ps.getPrimitives()
       self.failUnless(len(pp) == 8, "wrong number of items")
       self.failUnless(ps in pp, "ps missing in items.")
       self.failUnless(l0 in pp, "l0 missing in items.")
       self.failUnless(l1 in pp, "l1 missing in items.")
       self.failUnless(l2 in pp, "l2 missing in items.")
       self.failUnless(p0 in pp, "p0 missing in items.")
       self.failUnless(p1 in pp, "p1 missing in items.")
       self.failUnless(p3 in pp, "p3 missing in items.")
       self.failUnless(p4 in pp, "p4 missing in items.")
       
       ps.clearItems()
       self.failUnless(len(ps.getItems()) == 0, "cleaning items failed.")
          
   def test_PropertySet2D(self):
       p0=Point(0,0,0,0.1)
       p1=Point(10,0,0,0.2)
       p2=Point(10,10,0,0.3)
       p3=Point(0,10,3,0.4)
       p4=Point(5,5,0,0.001)
       p5=Point(7,5,0,0.001)
       p6=Point(5,7,0,0.001)

       l0=Line(p0,p1)
       l1=Line(p1,p2)
       l2=Line(p2,p3)
       l3=Line(p3,p0)

       l4=Line(p4,p5)
       l5=Line(p5,p6)
       l6=Line(p6,p4)

       cl=CurveLoop(l0,l1,l2,l3)
       h=CurveLoop(l4,l5,l6)

       s=PlaneSurface(cl,holes=[h])

       # create property set with dim:
       self.failUnlessRaises(TypeError,PropertySet,"test0",s, p0)

       ps=PropertySet("test0", s)
       self.failUnless(ps.getManifoldClass() == Manifold2D, "wrong manifold")
       self.failUnless(ps.getDim() == 2, "wrong dimension")
       self.failUnless(ps.getName() == "test0", "wrong name")
       ps.setName("test1")
       self.failUnless(ps.getName() == "test1", "wrong new name")

       self.failUnless(ps.getTag() == 19, "wrong tag")
       
       pp=ps.getPrimitives()
       self.failUnless(len(pp) == 18, "wrong number of items")
       self.failUnless(ps in pp, "ps missing in items.")
       self.failUnless(s in pp, "s missing in items.")
       self.failUnless(h in pp, "h missing in items.")
       self.failUnless(cl in pp, "cl missing in items.")
       self.failUnless(l0 in pp, "l0 missing in items.")
       self.failUnless(l1 in pp, "l1 missing in items.")
       self.failUnless(l2 in pp, "l2 missing in items.")
       self.failUnless(l3 in pp, "l3 missing in items.")
       self.failUnless(l4 in pp, "l4 missing in items.")
       self.failUnless(l5 in pp, "l5 missing in items.")
       self.failUnless(l6 in pp, "l6 missing in items.")
       self.failUnless(p0 in pp, "p0 missing in items.")
       self.failUnless(p1 in pp, "p1 missing in items.")
       self.failUnless(p3 in pp, "p3 missing in items.")
       self.failUnless(p4 in pp, "p4 missing in items.")
       self.failUnless(p5 in pp, "p5 missing in items.")
       self.failUnless(p6 in pp, "p6 missing in items.")

       self.failUnlessRaises(TypeError,ps.addItem, p0)
       pp=ps.getItems()
       self.failUnless(len(pp) == 1, "wrong number of items")
       self.failUnless(s in pp, "s missing in items.")
       
       ps.clearItems()
       self.failUnless(len(ps.getItems()) == 0, "cleaning items failed.")
       ps.addItem(s)
       self.failUnless(len(ps.getItems()) == 1, "missing added item")
       self.failUnless(s in ps.getItems(), "missing added item")

   def test_PropertySet3D(self):
       p0=Point(-2,-2,-2,0.1)
       p1=Point(2,-2,-2,0.1)
       p2=Point(-2,2,-2,0.1)
       p3=Point(2,2,-2,0.1)
       p4=Point(-2,-2,2,0.1)
       p5=Point(2,-2,2,0.1)
       p6=Point(-2,2,2,0.1)
       p7=Point(2,2,2,0.1)
       l01=Line(p0,p1)
       l15=Line(p1,p5)
       l54=Line(p5,p4)
       l40=Line(p4,p0)
       l23=Line(p2,p3)
       l37=Line(p3,p7)
       l76=Line(p7,p6)
       l62=Line(p6,p2)
       l13=Line(p1,p3)
       l57=Line(p5,p7)
       l02=Line(p0,p2)
       l46=Line(p4,p6)
       cl1=CurveLoop(l01,l15,l54,l40)
       s1=PlaneSurface(cl1)
       cl2=CurveLoop(l23,l37,l76,l62)
       s2=PlaneSurface(-cl2)
       cl3=CurveLoop(l13,l37,-l57,-l15)
       s3=PlaneSurface(cl3)
       cl4=CurveLoop(l46,l62,-l02,-l40)
       s4=PlaneSurface(-cl4)
       cl5=CurveLoop(-l01,l02,l23,-l13)
       s5=PlaneSurface(-cl5)
       cl6=CurveLoop(-l54,l57,l76,-l46)
       s6=PlaneSurface(-cl6)
       s_out=SurfaceLoop(s1,s2,s3,s4,s5,s6)

       p0_i=Point(-1,-1,-1,0.1)
       p1_i=Point(1,-1,-1,0.1)
       p2_i=Point(-1,1,-1,0.1)
       p3_i=Point(1,1,-1,0.1)
       p4_i=Point(-1,-1,1,0.1)
       p5_i=Point(1,-1,1,0.1)
       p6_i=Point(-1,1,1,0.1)
       p7_i=Point(1,1,1,0.1)
       l01_i=Line(p0_i,p1_i)
       l15_i=Line(p1_i,p5_i)
       l54_i=Line(p5_i,p4_i)
       l40_i=Line(p4_i,p0_i)
       l23_i=Line(p2_i,p3_i)
       l37_i=Line(p3_i,p7_i)
       l76_i=Line(p7_i,p6_i)
       l62_i=Line(p6_i,p2_i)
       l13_i=Line(p1_i,p3_i)
       l57_i=Line(p5_i,p7_i)
       l02_i=Line(p0_i,p2_i)
       l46_i=Line(p4_i,p6_i)
       cl1_i=CurveLoop(l01_i,l15_i,l54_i,l40_i)
       s1_i=PlaneSurface(cl1_i)
       cl2_i=CurveLoop(l23_i,l37_i,l76_i,l62_i)
       s2_i=PlaneSurface(-cl2_i)
       cl3_i=CurveLoop(l13_i,l37_i,-l57_i,-l15_i)
       s3_i=PlaneSurface(cl3_i)
       cl4_i=CurveLoop(l46_i,l62_i,-l02_i,-l40_i)
       s4_i=PlaneSurface(-cl4_i)
       cl5_i=CurveLoop(-l01_i,l02_i,l23_i,-l13_i)
       s5_i=PlaneSurface(-cl5_i)
       cl6_i=CurveLoop(-l54_i,l57_i,l76_i,-l46_i)
       s6_i=PlaneSurface(-cl6_i)
       s_inner=SurfaceLoop(s1_i,s2_i,s3_i,s4_i,s5_i,s6_i)

       v=Volume(s_out,holes=[s_inner])
          
       # create property set with dim:
       self.failUnlessRaises(TypeError,PropertySet,"test0",v, p0)

       ps=PropertySet("test0", v)
       self.failUnless(ps.getManifoldClass() == Manifold3D, "wrong manifold")
       self.failUnless(ps.getDim() == 3, "wrong dimension")
       self.failUnless(ps.getName() == "test0", "wrong name")
       ps.setName("test1")
       self.failUnless(ps.getName() == "test1", "wrong new name")

       self.failUnless(ps.getTag() == 69, "wrong tag")
       
       pp=ps.getPrimitives()
       self.failUnless(len(pp) == 68, "too many primitives.")
       self.failUnless(ps in pp, "ps is missing")
       self.failUnless(v in pp, "v is missing")
       self.failUnless(p0 in pp, "p0 is missing")
       self.failUnless(p1 in pp, "p1 is missing")
       self.failUnless(p2 in pp, "p2 is missing")
       self.failUnless(p3 in pp, "p3 is missing")
       self.failUnless(p4 in pp, "p4 is missing")
       self.failUnless(p5 in pp, "p5 is missing")
       self.failUnless(p6 in pp, "p6 is missing")
       self.failUnless(p7 in pp, "p7 is missing")
       self.failUnless(l01 in pp, "l01 is missing")
       self.failUnless(l15 in pp, "l15 is missing")
       self.failUnless(l54 in pp, "l54 is missing")
       self.failUnless(l40 in pp, "l40 is missing")
       self.failUnless(l23 in pp, "l23 is missing")
       self.failUnless(l37 in pp, "l37 is missing")
       self.failUnless(l76 in pp, "l76 is missing")
       self.failUnless(l62 in pp, "l62 is missing")
       self.failUnless(l13 in pp, "l13 is missing")
       self.failUnless(l57 in pp, "l57 is missing")
       self.failUnless(l02 in pp, "l02 is missing")
       self.failUnless(l46 in pp, "l46 is missing")
       self.failUnless(cl1 in pp, "cl1 is missing")
       self.failUnless(s1 in pp, "s1 is missing")
       self.failUnless(cl2 in pp, "cl2 is missing")
       self.failUnless(s2 in pp, "s2 is missing")
       self.failUnless(cl3 in pp, "cl3  is missing")
       self.failUnless(s3 in pp, "s3  is missing")
       self.failUnless(cl4 in pp, "cl4  is missing")
       self.failUnless(s4 in pp, "s4  is missing")
       self.failUnless(cl5 in pp, "cl5  is missing")
       self.failUnless(s5 in pp, "s5  is missing")
       self.failUnless(cl6 in pp, "cl6  is missing")
       self.failUnless(s6 in pp, "s6  is missing")
       self.failUnless(s_out in pp, "s_out is missing")
       self.failUnless(p0_i in pp, "p0_i is missing")
       self.failUnless(p1_i in pp, "p1_i is missing")
       self.failUnless(p2_i in pp, "p2_i is missing")
       self.failUnless(p3_i in pp, "p3_i is missing")
       self.failUnless(p4_i in pp, "p4_i is missing")
       self.failUnless(p5_i in pp, "p5_i is missing")
       self.failUnless(p6_i in pp, "p6_i is missing")
       self.failUnless(p7_i in pp, "p7_i is missing")
       self.failUnless(l01_i in pp, "l01_i is missing")
       self.failUnless(l15_i in pp, "l15_i is missing")
       self.failUnless(l54_i in pp, "l54_i is missing")
       self.failUnless(l40_i in pp, "l40_i is missing")
       self.failUnless(l23_i in pp, "l23_i is missing")
       self.failUnless(l37_i in pp, "l37_i is missing")
       self.failUnless(l76_i in pp, "l76_i is missing")
       self.failUnless(l62_i in pp, "l62_i is missing")
       self.failUnless(l13_i in pp, "l13_i is missing")
       self.failUnless(l57_i in pp, "l57_i is missing")
       self.failUnless(l02_i in pp, "l02_i is missing")
       self.failUnless(l46_i in pp, "l46_i is missing")
       self.failUnless(cl1_i in pp, "cl1_i is missing")
       self.failUnless(s1_i in pp, "s1_i is missing")
       self.failUnless(cl2_i in pp, "cl2_i is missing")
       self.failUnless(s2_i in pp, "s2_i is missing")
       self.failUnless(cl3_i in pp, "cl3_i  is missing")
       self.failUnless(s3_i in pp, "s3_i  is missing")
       self.failUnless(cl4_i in pp, "cl4_i  is missing")
       self.failUnless(s4_i in pp, "s4_i  is missing")
       self.failUnless(cl5_i in pp, "cl5_i  is missing")
       self.failUnless(s5_i in pp, "s5_i  is missing")
       self.failUnless(cl6_i in pp, "cl6_i  is missing")
       self.failUnless(s6_i in pp, "s6_i  is missing")
       self.failUnless(s_inner in pp, "s_inner  is missing")

       self.failUnlessRaises(TypeError,ps.addItem, p0)
       pp=ps.getItems()
       self.failUnless(len(pp) == 1, "wrong number of items")
       self.failUnless(v in pp, "s missing in items.")
       
       ps.clearItems()
       self.failUnless(len(ps.getItems()) == 0, "cleaning items failed.")
       ps.addItem(v)
       self.failUnless(len(ps.getItems()) == 1, "missing added item")
       self.failUnless(v in ps.getItems(), "missing added item")

class Test_PyCAD_Design(unittest.TestCase):
   def setUp(self):
         resetGlobalPrimitiveIdCounter()

   def test_tagMap(self):
       self.failUnlessRaises(TypeError,TagMap, { "x" : 5 } )
       self.failUnlessRaises(TypeError,TagMap, { 5 : 10 } )

       m=TagMap({ 5 : "x" })

       m.setMap(x=4,a=6)
       m.setMap(b=6,c=1)
       
       t=m.getTags()
       self.failUnless(len(t) == 4)
       self.failUnless(1 in t)
       self.failUnless(6 in t)
       self.failUnless(5 in t)
       self.failUnless(4 in t)
       self.failUnless(m.getTags("c") == [1])
       self.failUnless(m.getTags("x") == [4, 5])
       self.failUnless(m.getTags("b") == [6])

       self.failUnless(m.getName(1) == "c")
       self.failUnless(m.getName(4) == "x")
       self.failUnless(m.getName(5) == "x")
       self.failUnless(m.getName(6) == "b")

       t=m.getMapping()
       self.failUnless(len(t)==4)
       self.failUnless(t.has_key(1))
       self.failUnless(t.has_key(4))
       self.failUnless(t.has_key(5))
       self.failUnless(t.has_key(6))
       self.failUnless(t[1] == "c")
       self.failUnless(t[4] == "x")
       self.failUnless(t[5] == "x")
       self.failUnless(t[6] == "b")
       
       d=m.map(c=10, x = -10, b =60, default =1)
       self.failUnless(d[1] == 10)
       self.failUnless(d[4] == -10)
       self.failUnless(d[5] == -10)
       self.failUnless(d[6] == 60)

      
       d=m.map(c=10, x = -10, default =60)
       self.failUnless(d[1] == 10)
       self.failUnless(d[4] == -10)
       self.failUnless(d[5] == -10)
       self.failUnless(d[6] == 60)

       mxml=m.writeXML()
       m2=TagMap()
       m2.fillFromXML(mxml)
       t=m2.getMapping()
       self.failUnless(len(t)==4)
       self.failUnless(t.has_key(1))
       self.failUnless(t.has_key(4))
       self.failUnless(t.has_key(5))
       self.failUnless(t.has_key(6))
       self.failUnless(t[1] == "c")
       self.failUnless(t[4] == "x")
       self.failUnless(t[5] == "x")
       self.failUnless(t[6] == "b")

   def test_Design(self):
     
       d=Design0(dim=2, element_size=0.01, order=1, keep_files=False)
       # check dimension:
       self.failUnlessRaises(ValueError,d.setDim,4)
       d.setDim(3)
       self.failUnless(d.getDim() == 3)
       # check element order
       self.failUnlessRaises(ValueError,d.setElementOrder,4)
       d.setElementOrder(2)
       self.failUnless(d.getElementOrder() == 2)
       # check element size
       self.failUnlessRaises(ValueError,d.setElementSize,0)
       d.setElementSize(0.02)
       self.failUnless(d.getElementSize() == 0.02)
       # files:
       d.setKeepFilesOff()
       self.failUnless(not d.keepFiles())
       d.setKeepFilesOn()
       self.failUnless(d.keepFiles())
       # mesh handler:
       self.failUnlessRaises(NotImplementedError,d.getMeshHandler)

       p0=Point(0.,0.,0.)
       p1=Point(1.,0.,0.)
       p2=Point(1.,1.,0.)
       p3=Point(0.,1.,0.)
       l01=Line(p0,p1)
       l12=Line(p1,p2)
       l23=Line(p2,p3)
       l30=Line(p3,p0)
       c=CurveLoop(l01,l12,l23,l30)
       s=PlaneSurface(c)

       ps=PropertySet("XXXX",s)
       pl1=PropertySet("A",l01,l30)
       pl2=PropertySet("B",l12,l23)
       self.failUnlessRaises(TypeError,d.addItems,1.)
       d.addItems(ps, pl1, pl2)
       i=d.getItems()
       self.failUnless(isinstance(i,list))
       self.failUnless(len(i)==3)
       self.failUnless(ps in i)
       self.failUnless(pl1 in i)
       self.failUnless(pl2 in i)

       p=d.getAllPrimitives()
       self.failUnless(isinstance(p,list))
       self.failUnless(len(p)==13)
       self.failUnless(p0 in p)
       self.failUnless(p1 in p)
       self.failUnless(p2 in p)
       self.failUnless(p3 in p)
       self.failUnless(l01 in p)
       self.failUnless(l12 in p)
       self.failUnless(l23 in p)
       self.failUnless(l30 in p)
       self.failUnless(c in p)
       self.failUnless(s in p)
       self.failUnless(ps in p)
       self.failUnless(pl1 in p)
       self.failUnless(pl2 in p)
       # get tag maps:
       m1=d.getTagMap()
       mm1=m1.getMapping()
       self.failUnless(len(mm1)==3)
       self.failUnless(mm1[12]=="A")
       self.failUnless(mm1[13]=="B")
       self.failUnless(mm1[11]=="XXXX")
       # clear things:
       d.clearItems()
       i=d.getItems()
       self.failUnless(isinstance(i,list))
       self.failUnless(len(i)==0)


       
   def test_GMSH(self):
     
       d=GMSHDesign(dim=2, element_size=0.01, order=1, keep_files=False)

       script_name=d.getScriptFileName()
       self.failUnless(isinstance(script_name,str))
       self.failUnless(script_name.split(".")[-1] == "geo")
       script_name=PYCAD_WORKDIR+os.sep+"script.geo"
       d.setScriptFileName(script_name)
       self.failUnless(script_name == d.getScriptFileName())

       mesh_name=d.getMeshFileName()
       self.failUnless(isinstance(mesh_name,str))
       self.failUnless(mesh_name.split(".")[-1] == "msh")
       mesh_name=PYCAD_WORKDIR+os.sep+"mesh.msh"
       d.setMeshFileName(mesh_name)
       self.failUnless(mesh_name == d.getMeshFileName())
       
       d.setOptions(algorithm=d.TETGEN,optimize_quality=False,smoothing=4,curvature_based_element_size=True)
       cmd=d.getCommandString()
       self.failUnless("gmsh -2 -algo tetgen -clcurv -smooth 4 -v 0 -order 1 -o .%smesh.msh .%sscript.geo"%(os.sep,os.sep) == cmd)

       d.setOptions(optimize_quality=True,curvature_based_element_size=True)
       cmd=d.getCommandString()
       self.failUnless("gmsh -2 -algo iso -clcurv -smooth 1 -optimize -v 0 -order 1 -o .%smesh.msh .%sscript.geo"%(os.sep,os.sep) == cmd)

       p0=Point(0.,0.,0.)
       p1=Point(1.,0.,0.)
       p2=Point(1.,1.,0.)
       p3=Point(0.,1.,0.)
       l01=Line(p0,p1)
       l12=Line(p1,p2)
       l23=Line(p2,p3)
       l30=Line(p3,p0)
       c=CurveLoop(l01,l12,l23,l30)
       s=PlaneSurface(c)
       ps=PropertySet("XXXX",s)
       pl1=PropertySet("A",l01,l30)
       pl2=PropertySet("B",l12,l23)
       d.addItems(s,pl1,pl2)
       d.addItems(ps)
       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.01 };
Point(2) = {1.0 , 0.0, 0.0 , 0.01 };
Point(3) = {1.0 , 1.0, 0.0 , 0.01 };
Point(4) = {0.0 , 1.0, 0.0 , 0.01 };
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};
Physical Surface(11) = {10};
Physical Line(12) = {5, 8};
Physical Line(13) = {6, 7};
"""
       self.failUnless(scrpt == ref )

       
   def test_Triangle(self):
     
       d=TriangleDesign(dim=2, keep_files=False)

       script_name=d.getScriptFileName()
       self.failUnless(isinstance(script_name,str))
       self.failUnless(script_name.split(".")[-1] == "poly")
       script_name=PYCAD_WORKDIR+os.sep+"script.poly"
       d.setScriptFileName(script_name)
       self.failUnless(script_name == d.getScriptFileName())

       mesh_name=d.getMeshFileName()
       self.failUnless(isinstance(mesh_name,str))
       mesh_name=PYCAD_WORKDIR+os.sep+"mesh"
       d.setMeshFileName(mesh_name)
       self.failUnless(mesh_name == d.getMeshFileName())
       
       d.setOptions(cmdLineArgs="-Qpqa7.5")
       cmd=d.getCommandString()
       self.failUnless("triangle -Qpqa7.5 .%sscript.poly"%(os.sep) == cmd)

       p0=Point(0.,0.,0.)
       p1=Point(1.,0.,0.)
       p2=Point(1.,1.,0.)
       p3=Point(0.,1.,0.)
       l01=Line(p0,p1)
       l12=Line(p1,p2)
       l23=Line(p2,p3)
       l30=Line(p3,p0)
       c=CurveLoop(l01,l12,l23,l30)
       s=PlaneSurface(c)
       ps=PropertySet("XXXX",s)
       d.addItems(ps)

       scrpt=d.getScriptString()
       ref = \
"""# generated by esys.pycad
# vertices #
4 2 0 1
1 0.0 0.0 11
2 1.0 0.0 11
3 1.0 1.0 11
4 0.0 1.0 11
# segments #
4 1
1 1 2 11
2 2 3 11
3 3 4 11
4 4 1 11
# holes #
0
"""
       self.failUnless(scrpt == ref )

   def test_generate_Point(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       d.addItems(Point(1.,2.,3.,local_scale=9.))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
Point(1) = {1.0 , 2.0, 3.0 , 0.09 };
"""
       self.failUnless(scrpt == ref )
  
 
   def test_generate_Spline(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
 
        d.addItems(Spline(p0,p1,p2,p3))

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
Spline(5) = {1, 2, 3, 4};
"""
        self.failUnless(scrpt == ref )

   def test_generate_ReverseSpline(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
 
        CC0=Spline(p0,p1,p2,p3)
        d.addItems(-CC0)

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
Spline(5) = {1, 2, 3, 4};
"""
        self.failUnless(scrpt == ref )

   def test_generate_BezierCurve(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
 
        d.addItems(BezierCurve(p0,p1,p2,p3))

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
Bezier(5) = {1, 2, 3, 4};
"""
        self.failUnless(scrpt == ref )

   def test_generate_BSpline(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
 
        self.failUnlessRaises(ValueError,BSpline,p0)
        d.addItems(BSpline(p0,p1,p2,p3))

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
BSpline(5) = {1, 2, 3, 4};
"""
        self.failUnless(scrpt == ref )

   def test_generate_ReverseBSpline(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
 
        CC0=BSpline(p0,p1,p2,p3)
        d.addItems(-CC0)

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
BSpline(6) = {1, 2, 3, 4};
"""
        self.failUnless(scrpt == ref )

   def test_generate_LineSegment(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
 
        d.addItems(Line(p0,p1))

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Line(3) = {1, 2};
"""
        self.failUnless(scrpt == ref )

   def test_generate_ReverseLineSegment(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
 
        CC0=Line(p0,p1)
        d.addItems(-CC0)

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Line(3) = {1, 2};
"""
        self.failUnless(scrpt == ref )

   def test_generate_Arc(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        center=Point(0,0,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
 
        d.addItems(Arc(center,p_start,p_end))

        scrpt=d.getScriptString() 
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {1.0 , 2.0, 3.0 , 0.01 };
Circle(4) = {2, 1, 3};
"""
        self.failUnless(scrpt == ref )

   def test_generate_ReverseArc(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        center=Point(0,0,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
 
        CC0=Arc(center,p_start,p_end)
        d.addItems(-CC0)

        scrpt=d.getScriptString() 
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {1.0 , 2.0, 3.0 , 0.01 };
Circle(4) = {2, 1, 3};
"""
        self.failUnless(scrpt == ref )

   def test_generate_CurveLoop(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)

        l01=Line(p0,p1)
        l12=Arc(p3,p1,p2)
        l20=Spline(p2,p4,p0)

        lx=Line(p2,p3)
        ly=Line(p3,p1)

        d.addItems(CurveLoop(l01,l12,l20))

        scrpt=d.getScriptString() 
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
Point(5) = {1.0 , 2.0, 3.0 , 0.01 };
Line(6) = {1, 2};
Circle(7) = {2, 4, 3};
Spline(8) = {3, 5, 1};
Line Loop(11) = {6, 7, 8};
"""
        self.failUnless(scrpt == ref )

      
   def test_generate_ReverseCurveLoop(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)

        l01=Line(p0,p1)
        l12=Arc(p3,p1,p2)
        l20=Spline(p2,p4,p0)

        lx=Line(p2,p3)
        ly=Line(p3,p1)

        CC0=CurveLoop(l01,l20,l12)
        d.addItems(-CC0)

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
Point(5) = {1.0 , 2.0, 3.0 , 0.01 };
Line(6) = {1, 2};
Circle(7) = {2, 4, 3};
Spline(8) = {3, 5, 1};
Line Loop(11) = {6, 8, 7};
"""
        self.failUnless(scrpt == ref )

   def test_generate_RuledSurface(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        p5=Point(10,20,3)

        l01=Line(p0,p1)
        l12_1=Arc(p3,p1,p2)
        l12_2_1=Spline(p1,p3,p4)
        l12_2_2=Spline(p4,p5,p2)
        l12_3=Line(p1,p2)
        l20=Spline(p2,p4,p0)

        cl1=CurveLoop(l01,l12_1,l20) 
        cl2=CurveLoop(l01,l12_2_1,l12_2_2,l20)
        cl3=CurveLoop(l01,l12_3,l20)

        d.addItems(RuledSurface(cl1))

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
Point(5) = {1.0 , 2.0, 3.0 , 0.01 };
Line(7) = {1, 2};
Circle(8) = {2, 4, 3};
Spline(12) = {3, 5, 1};
Line Loop(13) = {7, 8, 12};
Ruled Surface(16) = {13};
"""
        self.failUnless(scrpt == ref )

   def test_generate_ReverseRuledSurface(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        p5=Point(10,20,3)

        l01=Line(p0,p1)
        l12_1=Arc(p3,p1,p2)
        l12_2_1=Spline(p1,p3,p4)
        l12_2_2=Spline(p4,p5,p2)
        l12_3=Line(p1,p2)
        l20=Spline(p2,p4,p0)

        cl1=CurveLoop(l01,l12_1,l20) 
        cl2=CurveLoop(l01,l12_2_1,l12_2_2,l20)
        cl3=CurveLoop(l01,l12_3,l20)

        CC0=RuledSurface(cl1)
        d.addItems(-CC0)

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {1.0 , 1.0, 1.0 , 0.002 };
Point(3) = {2.0 , 2.0, 2.0 , 0.003 };
Point(4) = {3.0 , 3.0, 3.0 , 0.004 };
Point(5) = {1.0 , 2.0, 3.0 , 0.01 };
Line(7) = {1, 2};
Circle(8) = {2, 4, 3};
Spline(12) = {3, 5, 1};
Line Loop(13) = {7, 8, 12};
Ruled Surface(16) = {13};
"""
        self.failUnless(scrpt == ref )

   def test_generate_PlaneSurface(self):
        d=GMSHDesign(dim=2, element_size=0.01)
        p0=Point(0,0,0,0.1)
        p1=Point(10,0,0,0.2)
        p2=Point(10,10,0,0.3)
        p3=Point(0,10,3,0.4)
        p4=Point(5,5,0,0.001)
        p5=Point(7,5,0,0.001)
        p6=Point(5,7,0,0.001)
        p7=Point(8,8,0,0.001)
        p8=Point(9,9,0,0.001)

        l0=Line(p0,p1)
        l1=Line(p1,p2)
        l2=Line(p2,p3)
        l3=Line(p3,p0)

        l9=Line(p1,p8)
        l10=Line(p8,p3)

        l4=Line(p4,p5)
        l5=Line(p5,p6)
        l6=Line(p6,p4)
        l7=Line(p6,p7)
        l8=Line(p7,p4)

        a1=Arc(p4,p3,p1)
        a2=Arc(p7,p5,p6)

        cl=CurveLoop(l0,l1,l2,l3)
        h=CurveLoop(l4,l5,l6)
        cl_s=CurveLoop(l0,l9,l10,l3)
        h2=CurveLoop(l4,l5,l7,l8)
        cl_a=CurveLoop(a1,l1,l2)
        h_a=CurveLoop(a2,l6,l4)

        d.addItems(PlaneSurface(cl,holes=[h]))

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {10.0 , 0.0, 0.0 , 0.002 };
Point(3) = {10.0 , 10.0, 0.0 , 0.003 };
Point(4) = {0.0 , 10.0, 3.0 , 0.004 };
Point(5) = {5.0 , 5.0, 0.0 , 1e-05 };
Point(6) = {7.0 , 5.0, 0.0 , 1e-05 };
Point(7) = {5.0 , 7.0, 0.0 , 1e-05 };
Line(10) = {1, 2};
Line(11) = {2, 3};
Line(12) = {3, 4};
Line(13) = {4, 1};
Line(16) = {5, 6};
Line(17) = {6, 7};
Line(18) = {7, 5};
Line Loop(23) = {10, 11, 12, 13};
Line Loop(24) = {16, 17, 18};
Plane Surface(29) = {23, 24};
"""
        self.failUnless(scrpt == ref )

   def test_generate_SurfaceLoop(self):
        d=GMSHDesign(dim=3, element_size=0.01)
        p0=Point( 0, 0, 0,0.1)
        p1=Point(10, 0, 0,0.1)
        p2=Point( 0,10, 0,0.1)
        p3=Point(10,10, 0,0.1)
        p4=Point( 0, 0,10,0.1)
        p5=Point(10, 0,10,0.1)
        p6=Point( 0,10,10,0.1)
        p7=Point(10,10,10,0.1)

        q0=Point( 4, 0, 4,0.1)
        q1=Point( 6, 0, 4,0.1)
        q2=Point( 4,10, 4,0.1)
        q3=Point( 6,10, 4,0.1)
        q4=Point( 4, 0, 6,0.1)
        q5=Point( 6, 0, 6,0.1)
        q6=Point( 4,10, 6,0.1)
        q7=Point( 6,10, 6,0.1)

        l01=Line(p0,p1)
        l15=Line(p1,p5)
        l54=Line(p5,p4)
        l04=Line(p4,p0)

        l13=Line(p1,p3)
        l37=Line(p3,p7)
        l75=Line(p7,p5)
        l67=Line(p6,p7)
        l26=Line(p2,p6)
        l32=Line(p3,p2)
        l20=Line(p2,p0)
        l46=Line(p4,p6)

        m01=Line(q0,q1)
        m15=Line(q1,q5)
        m54=Line(q5,q4)
        m40=Line(q4,q0)
        m23=Line(q2,q3)
        m37=Line(q3,q7)
        m76=Line(q7,q6)
        m62=Line(q6,q2)

        m02=Line(q0,q2)
        m13=Line(q1,q3)
        m46=Line(q4,q6)
        m57=Line(q5,q7)

        cl_l1=CurveLoop(l01,l15,l54,l04)
        cl_m1=CurveLoop(m01,m15,m54,m40)
        s1=PlaneSurface(cl_l1,holes=[cl_m1])
        s1_v=PlaneSurface(cl_l1)

        cl_l2=CurveLoop(-l15,l13,l37,l75)
        s2=PlaneSurface(cl_l2)

        cl_l3=CurveLoop(l32,-l37,l67,l26)
        cl_m3=CurveLoop(-m23,-m37,-m76,-m62)
        s3=PlaneSurface(cl_l3,holes=[cl_m3])
        s3_v=PlaneSurface(cl_l3)
     
     
        cl_l4=CurveLoop(l20,-l26,l46,-l04)
        s4=PlaneSurface(cl_l4)

        cl_l5=CurveLoop(l32,l20,l01,l13)
        s5=PlaneSurface(-cl_l5)

        cl_l6=CurveLoop(l67,l75,l54,l46)
        s6=PlaneSurface(-cl_l6)

        cl_m7=CurveLoop(m13,m37,-m57,-m15)
        s7=PlaneSurface(cl_m7)
        
        cl_m8=CurveLoop(m57,m76,-m46,-m54)
        s8=PlaneSurface(cl_m8)

        cl_m9=CurveLoop(m46,m62,-m02,-m40)
        s9=PlaneSurface(cl_m9)

        cl_m10=CurveLoop(-m01,m02,m23,-m13)
        s10=PlaneSurface(cl_m10)

        d.addItems(SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10))

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {10.0 , 0.0, 0.0 , 0.001 };
Point(3) = {0.0 , 10.0, 0.0 , 0.001 };
Point(4) = {10.0 , 10.0, 0.0 , 0.001 };
Point(5) = {0.0 , 0.0, 10.0 , 0.001 };
Point(6) = {10.0 , 0.0, 10.0 , 0.001 };
Point(7) = {0.0 , 10.0, 10.0 , 0.001 };
Point(8) = {10.0 , 10.0, 10.0 , 0.001 };
Point(9) = {4.0 , 0.0, 4.0 , 0.001 };
Point(10) = {6.0 , 0.0, 4.0 , 0.001 };
Point(11) = {4.0 , 10.0, 4.0 , 0.001 };
Point(12) = {6.0 , 10.0, 4.0 , 0.001 };
Point(13) = {4.0 , 0.0, 6.0 , 0.001 };
Point(14) = {6.0 , 0.0, 6.0 , 0.001 };
Point(15) = {4.0 , 10.0, 6.0 , 0.001 };
Point(16) = {6.0 , 10.0, 6.0 , 0.001 };
Line(17) = {1, 2};
Line(18) = {2, 6};
Line(19) = {6, 5};
Line(20) = {5, 1};
Line(21) = {2, 4};
Line(22) = {4, 8};
Line(23) = {8, 6};
Line(24) = {7, 8};
Line(25) = {3, 7};
Line(26) = {4, 3};
Line(27) = {3, 1};
Line(28) = {5, 7};
Line(29) = {9, 10};
Line(30) = {10, 14};
Line(31) = {14, 13};
Line(32) = {13, 9};
Line(33) = {11, 12};
Line(34) = {12, 16};
Line(35) = {16, 15};
Line(36) = {15, 11};
Line(37) = {9, 11};
Line(38) = {10, 12};
Line(39) = {13, 15};
Line(40) = {14, 16};
Line Loop(41) = {17, 18, 19, 20};
Line Loop(42) = {29, 30, 31, 32};
Plane Surface(43) = {41, 42};
Line Loop(45) = {-18, 21, 22, 23};
Plane Surface(46) = {45};
Line Loop(47) = {26, -22, 24, 25};
Line Loop(48) = {-33, -34, -35, -36};
Plane Surface(49) = {47, 48};
Line Loop(51) = {27, -25, 28, -20};
Plane Surface(52) = {51};
Line Loop(53) = {26, 27, 17, 21};
Plane Surface(54) = {-53};
Line Loop(55) = {24, 23, 19, 28};
Plane Surface(56) = {-55};
Line Loop(57) = {38, 34, -40, -30};
Plane Surface(58) = {57};
Line Loop(59) = {40, 35, -39, -31};
Plane Surface(60) = {59};
Line Loop(61) = {39, 36, -37, -32};
Plane Surface(62) = {61};
Line Loop(63) = {-29, 37, 33, -38};
Plane Surface(64) = {63};
Surface Loop(65) = {43, 46, 49, 52, 54, 56, -58, -60, -62, -64};
"""
        self.failUnless(scrpt == ref )

   def test_generate_ReverseSurfaceLoop(self):
        d=GMSHDesign(dim=3, element_size=0.01)
        p0=Point( 0, 0, 0,0.1)
        p1=Point(10, 0, 0,0.1)
        p2=Point( 0,10, 0,0.1)
        p3=Point(10,10, 0,0.1)
        p4=Point( 0, 0,10,0.1)
        p5=Point(10, 0,10,0.1)
        p6=Point( 0,10,10,0.1)
        p7=Point(10,10,10,0.1)

        q0=Point( 4, 0, 4,0.1)
        q1=Point( 6, 0, 4,0.1)
        q2=Point( 4,10, 4,0.1)
        q3=Point( 6,10, 4,0.1)
        q4=Point( 4, 0, 6,0.1)
        q5=Point( 6, 0, 6,0.1)
        q6=Point( 4,10, 6,0.1)
        q7=Point( 6,10, 6,0.1)

        l01=Line(p0,p1)
        l15=Line(p1,p5)
        l54=Line(p5,p4)
        l04=Line(p4,p0)

        l13=Line(p1,p3)
        l37=Line(p3,p7)
        l75=Line(p7,p5)
        l67=Line(p6,p7)
        l26=Line(p2,p6)
        l32=Line(p3,p2)
        l20=Line(p2,p0)
        l46=Line(p4,p6)

        m01=Line(q0,q1)
        m15=Line(q1,q5)
        m54=Line(q5,q4)
        m40=Line(q4,q0)
        m23=Line(q2,q3)
        m37=Line(q3,q7)
        m76=Line(q7,q6)
        m62=Line(q6,q2)

        m02=Line(q0,q2)
        m13=Line(q1,q3)
        m46=Line(q4,q6)
        m57=Line(q5,q7)

        cl_l1=CurveLoop(l01,l15,l54,l04)
        cl_m1=CurveLoop(m01,m15,m54,m40)
        s1=PlaneSurface(cl_l1,holes=[cl_m1])
        s1_v=PlaneSurface(cl_l1)

        cl_l2=CurveLoop(-l15,l13,l37,l75)
        s2=PlaneSurface(cl_l2)

        cl_l3=CurveLoop(l32,-l37,l67,l26)
        cl_m3=CurveLoop(-m23,-m37,-m76,-m62)
        s3=PlaneSurface(cl_l3,holes=[cl_m3])
        s3_v=PlaneSurface(cl_l3)
     
     
        cl_l4=CurveLoop(l20,-l26,l46,-l04)
        s4=PlaneSurface(cl_l4)

        cl_l5=CurveLoop(l32,l20,l01,l13)
        s5=PlaneSurface(-cl_l5)

        cl_l6=CurveLoop(l67,l75,l54,l46)
        s6=PlaneSurface(-cl_l6)

        cl_m7=CurveLoop(m13,m37,-m57,-m15)
        s7=PlaneSurface(cl_m7)
        
        cl_m8=CurveLoop(m57,m76,-m46,-m54)
        s8=PlaneSurface(cl_m8)

        cl_m9=CurveLoop(m46,m62,-m02,-m40)
        s9=PlaneSurface(cl_m9)

        cl_m10=CurveLoop(-m01,m02,m23,-m13)
        s10=PlaneSurface(cl_m10)

        CC0=SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)
        s=-CC0
        d.addItems(-CC0)

        scrpt=d.getScriptString()
        ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {10.0 , 0.0, 0.0 , 0.001 };
Point(3) = {0.0 , 10.0, 0.0 , 0.001 };
Point(4) = {10.0 , 10.0, 0.0 , 0.001 };
Point(5) = {0.0 , 0.0, 10.0 , 0.001 };
Point(6) = {10.0 , 0.0, 10.0 , 0.001 };
Point(7) = {0.0 , 10.0, 10.0 , 0.001 };
Point(8) = {10.0 , 10.0, 10.0 , 0.001 };
Point(9) = {4.0 , 0.0, 4.0 , 0.001 };
Point(10) = {6.0 , 0.0, 4.0 , 0.001 };
Point(11) = {4.0 , 10.0, 4.0 , 0.001 };
Point(12) = {6.0 , 10.0, 4.0 , 0.001 };
Point(13) = {4.0 , 0.0, 6.0 , 0.001 };
Point(14) = {6.0 , 0.0, 6.0 , 0.001 };
Point(15) = {4.0 , 10.0, 6.0 , 0.001 };
Point(16) = {6.0 , 10.0, 6.0 , 0.001 };
Line(17) = {1, 2};
Line(18) = {2, 6};
Line(19) = {6, 5};
Line(20) = {5, 1};
Line(21) = {2, 4};
Line(22) = {4, 8};
Line(23) = {8, 6};
Line(24) = {7, 8};
Line(25) = {3, 7};
Line(26) = {4, 3};
Line(27) = {3, 1};
Line(28) = {5, 7};
Line(29) = {9, 10};
Line(30) = {10, 14};
Line(31) = {14, 13};
Line(32) = {13, 9};
Line(33) = {11, 12};
Line(34) = {12, 16};
Line(35) = {16, 15};
Line(36) = {15, 11};
Line(37) = {9, 11};
Line(38) = {10, 12};
Line(39) = {13, 15};
Line(40) = {14, 16};
Line Loop(41) = {17, 18, 19, 20};
Line Loop(42) = {29, 30, 31, 32};
Plane Surface(43) = {41, 42};
Line Loop(45) = {-18, 21, 22, 23};
Plane Surface(46) = {45};
Line Loop(47) = {26, -22, 24, 25};
Line Loop(48) = {-33, -34, -35, -36};
Plane Surface(49) = {47, 48};
Line Loop(51) = {27, -25, 28, -20};
Plane Surface(52) = {51};
Line Loop(53) = {26, 27, 17, 21};
Plane Surface(54) = {-53};
Line Loop(55) = {24, 23, 19, 28};
Plane Surface(56) = {-55};
Line Loop(57) = {38, 34, -40, -30};
Plane Surface(58) = {57};
Line Loop(59) = {40, 35, -39, -31};
Plane Surface(60) = {59};
Line Loop(61) = {39, 36, -37, -32};
Plane Surface(62) = {61};
Line Loop(63) = {-29, 37, 33, -38};
Plane Surface(64) = {63};
Surface Loop(65) = {43, 46, 49, 52, 54, 56, -58, -60, -62, -64};
"""
        self.failUnless(scrpt == ref )

   def test_generate_Volume(self):
        d=GMSHDesign(dim=3, element_size=0.01)
        p0=Point(-2,-2,-2,0.1)
        p1=Point(2,-2,-2,0.1)
        p2=Point(-2,2,-2,0.1)
        p3=Point(2,2,-2,0.1)
        p4=Point(-2,-2,2,0.1)
        p5=Point(2,-2,2,0.1)
        p6=Point(-2,2,2,0.1)
        p7=Point(2,2,2,0.1)
        l01=Line(p0,p1)
        l15=Line(p1,p5)
        l54=Line(p5,p4)
        l40=Line(p4,p0)
        l23=Line(p2,p3)
        l37=Line(p3,p7)
        l76=Line(p7,p6)
        l62=Line(p6,p2)
        l13=Line(p1,p3)
        l57=Line(p5,p7)
        l02=Line(p0,p2)
        l46=Line(p4,p6)
        cl1=CurveLoop(l01,l15,l54,l40)
        s1=PlaneSurface(cl1)
        cl2=CurveLoop(l23,l37,l76,l62)
        s2=PlaneSurface(-cl2)
        cl3=CurveLoop(l13,l37,-l57,-l15)
        s3=PlaneSurface(cl3)
        cl4=CurveLoop(l46,l62,-l02,-l40)
        s4=PlaneSurface(-cl4)
        cl5=CurveLoop(-l01,l02,l23,-l13)
        s5=PlaneSurface(-cl5)
        cl6=CurveLoop(-l54,l57,l76,-l46)
        s6=PlaneSurface(-cl6)
        s_out=SurfaceLoop(s1,s2,s3,s4,s5,s6)
 
        p0_i=Point(-1,-1,-1,0.1)
        p1_i=Point(1,-1,-1,0.1)
        p2_i=Point(-1,1,-1,0.1)
        p3_i=Point(1,1,-1,0.1)
        p4_i=Point(-1,-1,1,0.1)
        p5_i=Point(1,-1,1,0.1)
        p6_i=Point(-1,1,1,0.1)
        p7_i=Point(1,1,1,0.1)
        l01_i=Line(p0_i,p1_i)
        l15_i=Line(p1_i,p5_i)
        l54_i=Line(p5_i,p4_i)
        l40_i=Line(p4_i,p0_i)
        l23_i=Line(p2_i,p3_i)
        l37_i=Line(p3_i,p7_i)
        l76_i=Line(p7_i,p6_i)
        l62_i=Line(p6_i,p2_i)
        l13_i=Line(p1_i,p3_i)
        l57_i=Line(p5_i,p7_i)
        l02_i=Line(p0_i,p2_i)
        l46_i=Line(p4_i,p6_i)
        cl1_i=CurveLoop(l01_i,l15_i,l54_i,l40_i)
        s1_i=PlaneSurface(cl1_i)
        cl2_i=CurveLoop(l23_i,l37_i,l76_i,l62_i)
        s2_i=PlaneSurface(-cl2_i)
        cl3_i=CurveLoop(l13_i,l37_i,-l57_i,-l15_i)
        s3_i=PlaneSurface(cl3_i)
        cl4_i=CurveLoop(l46_i,l62_i,-l02_i,-l40_i)
        s4_i=PlaneSurface(-cl4_i)
        cl5_i=CurveLoop(-l01_i,l02_i,l23_i,-l13_i)
        s5_i=PlaneSurface(-cl5_i)
        cl6_i=CurveLoop(-l54_i,l57_i,l76_i,-l46_i)
        s6_i=PlaneSurface(-cl6_i)
        s_inner=SurfaceLoop(s1_i,s2_i,s3_i,s4_i,s5_i,s6_i)

        d.addItems(Volume(s_out,holes=[s_inner]))

        scrpt=d.getScriptString(); 
        ref = \
"""// generated by esys.pycad
Point(1) = {-2.0 , -2.0, -2.0 , 0.001 };
Point(2) = {2.0 , -2.0, -2.0 , 0.001 };
Point(3) = {-2.0 , 2.0, -2.0 , 0.001 };
Point(4) = {2.0 , 2.0, -2.0 , 0.001 };
Point(5) = {-2.0 , -2.0, 2.0 , 0.001 };
Point(6) = {2.0 , -2.0, 2.0 , 0.001 };
Point(7) = {-2.0 , 2.0, 2.0 , 0.001 };
Point(8) = {2.0 , 2.0, 2.0 , 0.001 };
Line(9) = {1, 2};
Line(10) = {2, 6};
Line(11) = {6, 5};
Line(12) = {5, 1};
Line(13) = {3, 4};
Line(14) = {4, 8};
Line(15) = {8, 7};
Line(16) = {7, 3};
Line(17) = {2, 4};
Line(18) = {6, 8};
Line(19) = {1, 3};
Line(20) = {5, 7};
Line Loop(21) = {9, 10, 11, 12};
Plane Surface(22) = {21};
Line Loop(23) = {13, 14, 15, 16};
Plane Surface(24) = {-23};
Line Loop(25) = {17, 14, -18, -10};
Plane Surface(26) = {25};
Line Loop(27) = {20, 16, -19, -12};
Plane Surface(28) = {-27};
Line Loop(29) = {-9, 19, 13, -17};
Plane Surface(30) = {-29};
Line Loop(31) = {-11, 18, 15, -20};
Plane Surface(32) = {-31};
Surface Loop(33) = {22, 24, 26, 28, 30, 32};
Point(34) = {-1.0 , -1.0, -1.0 , 0.001 };
Point(35) = {1.0 , -1.0, -1.0 , 0.001 };
Point(36) = {-1.0 , 1.0, -1.0 , 0.001 };
Point(37) = {1.0 , 1.0, -1.0 , 0.001 };
Point(38) = {-1.0 , -1.0, 1.0 , 0.001 };
Point(39) = {1.0 , -1.0, 1.0 , 0.001 };
Point(40) = {-1.0 , 1.0, 1.0 , 0.001 };
Point(41) = {1.0 , 1.0, 1.0 , 0.001 };
Line(42) = {34, 35};
Line(43) = {35, 39};
Line(44) = {39, 38};
Line(45) = {38, 34};
Line(46) = {36, 37};
Line(47) = {37, 41};
Line(48) = {41, 40};
Line(49) = {40, 36};
Line(50) = {35, 37};
Line(51) = {39, 41};
Line(52) = {34, 36};
Line(53) = {38, 40};
Line Loop(54) = {42, 43, 44, 45};
Plane Surface(55) = {54};
Line Loop(56) = {46, 47, 48, 49};
Plane Surface(57) = {-56};
Line Loop(58) = {50, 47, -51, -43};
Plane Surface(59) = {58};
Line Loop(60) = {53, 49, -52, -45};
Plane Surface(61) = {-60};
Line Loop(62) = {-42, 52, 46, -50};
Plane Surface(63) = {-62};
Line Loop(64) = {-44, 51, 48, -53};
Plane Surface(65) = {-64};
Surface Loop(66) = {55, 57, 59, 61, 63, 65};
Volume(67) = {33, 66};
"""
        self.failUnless(scrpt == ref )

   def test_generate_PropertySet0D(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       p0=Point(1.,2.,3.,local_scale=9.)
       p1=Point(0.,0.,0.,local_scale=9.)
       p3=Point(8.,6,6,local_scale=9.)
       p4=Point(8.,6,-6,local_scale=9.)
       
       d.addItems(PropertySet("test0",p0, p1, p3, p4))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
Point(1) = {1.0 , 2.0, 3.0 , 0.09 };
Point(2) = {0.0 , 0.0, 0.0 , 0.09 };
Point(3) = {8.0 , 6.0, 6.0 , 0.09 };
Point(4) = {8.0 , 6.0, -6.0 , 0.09 };
Physical Point(5) = {1, 2, 3, 4};
"""
       self.failUnless(scrpt == ref )
       
   def test_generate_PropertySet1D(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       p0=Point(1.,2.,3.,local_scale=9.)
       p1=Point(0.,0.,0.,local_scale=9.)
       p3=Point(8.,6,6,local_scale=9.)
       p4=Point(8.,6,-6,local_scale=9.)
       
       l0=Line(p0,p1)
       l1=Arc(p3,p1,p4)
       l2=Line(p4,p0)

       d.addItems(PropertySet("test0", l0, l1, l2))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
Point(1) = {1.0 , 2.0, 3.0 , 0.09 };
Point(2) = {0.0 , 0.0, 0.0 , 0.09 };
Point(3) = {8.0 , 6.0, 6.0 , 0.09 };
Point(4) = {8.0 , 6.0, -6.0 , 0.09 };
Line(5) = {1, 2};
Circle(6) = {2, 3, 4};
Line(7) = {4, 1};
Physical Line(8) = {5, 6, 7};
"""
       self.failUnless(scrpt == ref )
 
   def test_generate_PropertySet2D(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       p0=Point(0,0,0,0.1)
       p1=Point(10,0,0,0.2)
       p2=Point(10,10,0,0.3)
       p3=Point(0,10,3,0.4)
       p4=Point(5,5,0,0.001)
       p5=Point(7,5,0,0.001)
       p6=Point(5,7,0,0.001)

       l0=Line(p0,p1)
       l1=Line(p1,p2)
       l2=Line(p2,p3)
       l3=Line(p3,p0)

       l4=Line(p4,p5)
       l5=Line(p5,p6)
       l6=Line(p6,p4)

       cl=CurveLoop(l0,l1,l2,l3)
       h=CurveLoop(l4,l5,l6)

       s=PlaneSurface(cl,holes=[h])

       d.addItems(PropertySet("test0", s))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
Point(1) = {0.0 , 0.0, 0.0 , 0.001 };
Point(2) = {10.0 , 0.0, 0.0 , 0.002 };
Point(3) = {10.0 , 10.0, 0.0 , 0.003 };
Point(4) = {0.0 , 10.0, 3.0 , 0.004 };
Point(5) = {5.0 , 5.0, 0.0 , 1e-05 };
Point(6) = {7.0 , 5.0, 0.0 , 1e-05 };
Point(7) = {5.0 , 7.0, 0.0 , 1e-05 };
Line(8) = {1, 2};
Line(9) = {2, 3};
Line(10) = {3, 4};
Line(11) = {4, 1};
Line(12) = {5, 6};
Line(13) = {6, 7};
Line(14) = {7, 5};
Line Loop(15) = {8, 9, 10, 11};
Line Loop(16) = {12, 13, 14};
Plane Surface(17) = {15, 16};
Physical Surface(18) = {17};
"""
       self.failUnless(scrpt == ref )

   def test_generate_PropertySet3D(self):
       d=GMSHDesign(dim=3, element_size=0.01)
       p0=Point(-2,-2,-2,0.1)
       p1=Point(2,-2,-2,0.1)
       p2=Point(-2,2,-2,0.1)
       p3=Point(2,2,-2,0.1)
       p4=Point(-2,-2,2,0.1)
       p5=Point(2,-2,2,0.1)
       p6=Point(-2,2,2,0.1)
       p7=Point(2,2,2,0.1)
       l01=Line(p0,p1)
       l15=Line(p1,p5)
       l54=Line(p5,p4)
       l40=Line(p4,p0)
       l23=Line(p2,p3)
       l37=Line(p3,p7)
       l76=Line(p7,p6)
       l62=Line(p6,p2)
       l13=Line(p1,p3)
       l57=Line(p5,p7)
       l02=Line(p0,p2)
       l46=Line(p4,p6)
       cl1=CurveLoop(l01,l15,l54,l40)
       s1=PlaneSurface(cl1)
       cl2=CurveLoop(l23,l37,l76,l62)
       s2=PlaneSurface(-cl2)
       cl3=CurveLoop(l13,l37,-l57,-l15)
       s3=PlaneSurface(cl3)
       cl4=CurveLoop(l46,l62,-l02,-l40)
       s4=PlaneSurface(-cl4)
       cl5=CurveLoop(-l01,l02,l23,-l13)
       s5=PlaneSurface(-cl5)
       cl6=CurveLoop(-l54,l57,l76,-l46)
       s6=PlaneSurface(-cl6)
       s_out=SurfaceLoop(s1,s2,s3,s4,s5,s6)

       p0_i=Point(-1,-1,-1,0.1)
       p1_i=Point(1,-1,-1,0.1)
       p2_i=Point(-1,1,-1,0.1)
       p3_i=Point(1,1,-1,0.1)
       p4_i=Point(-1,-1,1,0.1)
       p5_i=Point(1,-1,1,0.1)
       p6_i=Point(-1,1,1,0.1)
       p7_i=Point(1,1,1,0.1)
       l01_i=Line(p0_i,p1_i)
       l15_i=Line(p1_i,p5_i)
       l54_i=Line(p5_i,p4_i)
       l40_i=Line(p4_i,p0_i)
       l23_i=Line(p2_i,p3_i)
       l37_i=Line(p3_i,p7_i)
       l76_i=Line(p7_i,p6_i)
       l62_i=Line(p6_i,p2_i)
       l13_i=Line(p1_i,p3_i)
       l57_i=Line(p5_i,p7_i)
       l02_i=Line(p0_i,p2_i)
       l46_i=Line(p4_i,p6_i)
       cl1_i=CurveLoop(l01_i,l15_i,l54_i,l40_i)
       s1_i=PlaneSurface(cl1_i)
       cl2_i=CurveLoop(l23_i,l37_i,l76_i,l62_i)
       s2_i=PlaneSurface(-cl2_i)
       cl3_i=CurveLoop(l13_i,l37_i,-l57_i,-l15_i)
       s3_i=PlaneSurface(cl3_i)
       cl4_i=CurveLoop(l46_i,l62_i,-l02_i,-l40_i)
       s4_i=PlaneSurface(-cl4_i)
       cl5_i=CurveLoop(-l01_i,l02_i,l23_i,-l13_i)
       s5_i=PlaneSurface(-cl5_i)
       cl6_i=CurveLoop(-l54_i,l57_i,l76_i,-l46_i)
       s6_i=PlaneSurface(-cl6_i)
       s_inner=SurfaceLoop(s1_i,s2_i,s3_i,s4_i,s5_i,s6_i)

       v=Volume(s_out,holes=[s_inner])
       d.addItems(PropertySet("test0", v))
       
       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
Point(1) = {-2.0 , -2.0, -2.0 , 0.001 };
Point(2) = {2.0 , -2.0, -2.0 , 0.001 };
Point(3) = {-2.0 , 2.0, -2.0 , 0.001 };
Point(4) = {2.0 , 2.0, -2.0 , 0.001 };
Point(5) = {-2.0 , -2.0, 2.0 , 0.001 };
Point(6) = {2.0 , -2.0, 2.0 , 0.001 };
Point(7) = {-2.0 , 2.0, 2.0 , 0.001 };
Point(8) = {2.0 , 2.0, 2.0 , 0.001 };
Line(9) = {1, 2};
Line(10) = {2, 6};
Line(11) = {6, 5};
Line(12) = {5, 1};
Line(13) = {3, 4};
Line(14) = {4, 8};
Line(15) = {8, 7};
Line(16) = {7, 3};
Line(17) = {2, 4};
Line(18) = {6, 8};
Line(19) = {1, 3};
Line(20) = {5, 7};
Line Loop(21) = {9, 10, 11, 12};
Plane Surface(22) = {21};
Line Loop(23) = {13, 14, 15, 16};
Plane Surface(24) = {-23};
Line Loop(25) = {17, 14, -18, -10};
Plane Surface(26) = {25};
Line Loop(27) = {20, 16, -19, -12};
Plane Surface(28) = {-27};
Line Loop(29) = {-9, 19, 13, -17};
Plane Surface(30) = {-29};
Line Loop(31) = {-11, 18, 15, -20};
Plane Surface(32) = {-31};
Surface Loop(33) = {22, 24, 26, 28, 30, 32};
Point(34) = {-1.0 , -1.0, -1.0 , 0.001 };
Point(35) = {1.0 , -1.0, -1.0 , 0.001 };
Point(36) = {-1.0 , 1.0, -1.0 , 0.001 };
Point(37) = {1.0 , 1.0, -1.0 , 0.001 };
Point(38) = {-1.0 , -1.0, 1.0 , 0.001 };
Point(39) = {1.0 , -1.0, 1.0 , 0.001 };
Point(40) = {-1.0 , 1.0, 1.0 , 0.001 };
Point(41) = {1.0 , 1.0, 1.0 , 0.001 };
Line(42) = {34, 35};
Line(43) = {35, 39};
Line(44) = {39, 38};
Line(45) = {38, 34};
Line(46) = {36, 37};
Line(47) = {37, 41};
Line(48) = {41, 40};
Line(49) = {40, 36};
Line(50) = {35, 37};
Line(51) = {39, 41};
Line(52) = {34, 36};
Line(53) = {38, 40};
Line Loop(54) = {42, 43, 44, 45};
Plane Surface(55) = {54};
Line Loop(56) = {46, 47, 48, 49};
Plane Surface(57) = {-56};
Line Loop(58) = {50, 47, -51, -43};
Plane Surface(59) = {58};
Line Loop(60) = {53, 49, -52, -45};
Plane Surface(61) = {-60};
Line Loop(62) = {-42, 52, 46, -50};
Plane Surface(63) = {-62};
Line Loop(64) = {-44, 51, 48, -53};
Plane Surface(65) = {-64};
Surface Loop(66) = {55, 57, 59, 61, 63, 65};
Volume(67) = {33, 66};
Physical Volume(68) = {67};
"""
       self.failUnless(scrpt == ref )

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_PyCAD_Transformations))
   suite.addTest(unittest.makeSuite(Test_PyCAD_Primitives))
   suite.addTest(unittest.makeSuite(Test_PyCAD_Design))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
