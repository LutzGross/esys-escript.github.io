# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import os
import sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import math
import numpy
from esys.pycad import *
from esys.pycad.design import AbstractDesign
from esys.pycad.gmsh import Design as GMSHDesign
from esys.pycad.extras import layer_cake

try:
     PYCAD_TEST_DATA=os.environ['PYCAD_TEST_DATA']
except KeyError:
     PYCAD_TEST_DATA='.'

try:
     PYCAD_WORKDIR=os.environ['PYCAD_WORKDIR']
except KeyError:
     PYCAD_WORKDIR='.'

#PYCAD_TEST_MESH_PATH=PYCAD_TEST_DATA+os.sep+"data_meshes"+os.sep
#PYCAD_WORKDIR_PATH=PYCAD_WORKDIR+os.sep

def _cross(x, y):
    return numpy.array([x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]])


class Test_PyCAD_Transformations(unittest.TestCase):
   ABS_TOL=1.e-8
   def __distance(self,x,y):
       return math.sqrt(numpy.dot(x-y,x-y))
   def test_Translation_x(self):
        t=Translation([1,0,0])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([2,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([1,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([1,0,1]))<self.ABS_TOL,"s2 is wrong.")
   def test_Translation_y(self):
        t=Translation([0,1,0])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1,1,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,2,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,1,1]))<self.ABS_TOL,"s2 is wrong.")
   def test_Translation_z(self):
        t=Translation([0,0,1])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1,0,1]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,1,1]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_0_two(self):
        t=Dilation(2.)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([2,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,2,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_0_half(self):
        t=Dilation(0.5)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([0.5,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,0.5,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,0.5]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_x_two(self):
        t=Dilation(2.,[1.,0.,0.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s0_1=t([0,0,0])
        self.assertTrue(isinstance(s0_1,numpy.ndarray),"s0_1 is not an ndarray object.")
        self.assertTrue(self.__distance(s0_1,numpy.array([-1.,0,0]))<self.ABS_TOL,"s0_1 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([-1,2,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([-1.,0,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_x_half(self):
        t=Dilation(0.5,[1.,0.,0.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0]))<self.ABS_TOL,"s0 is wrong.")
        s0_1=t([0,0,0])
        self.assertTrue(isinstance(s0_1,numpy.ndarray),"s0_1 is not an ndarray object.")
        self.assertTrue(self.__distance(s0_1,numpy.array([.5,0,0]))<self.ABS_TOL,"s0_1 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.5,0.5,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.5,0,0.5]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_y_two(self):
        t=Dilation(2.,[0.,1.,0.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([2.,-1.,0]))<self.ABS_TOL,"s0 is wrong.")
        s1_1=t([0,0,0])
        self.assertTrue(isinstance(s1_1,numpy.ndarray),"s1_1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1_1,numpy.array([0.,-1.,0]))<self.ABS_TOL,"s1_1 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,1.,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,-1.,2]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_y_half(self):
        t=Dilation(0.5,[0.,1.,0.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([0.5,0.5,0]))<self.ABS_TOL,"s0 is wrong.")
        s1_1=t([0,0,0])
        self.assertTrue(isinstance(s1_1,numpy.ndarray),"s1_1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1_1,numpy.array([0,0.5,0]))<self.ABS_TOL,"s1_1 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,1.,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0.5,0.5]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_z_two(self):
        t=Dilation(2.,[0.,0.,1.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([2.,0.,-1.]))<self.ABS_TOL,"s0 is wrong.")
        s2_1=t([0,0,0])
        self.assertTrue(isinstance(s2_1,numpy.ndarray),"s2_1 is not an ndarray object.")
        self.assertTrue(self.__distance(s2_1,numpy.array([0.,0.,-1.]))<self.ABS_TOL,"s2_1 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,2.,-1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0.,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Dilation_z_half(self):
        t=Dilation(0.5,[0.,0.,1.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([0.5,0.,0.5]))<self.ABS_TOL,"s0 is wrong.")
        s2_1=t([0,0,0])
        self.assertTrue(isinstance(s2_1,numpy.ndarray),"s2_1 is not an ndarray object.")
        self.assertTrue(self.__distance(s2_1,numpy.array([0,0,0.5]))<self.ABS_TOL,"s2_1 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,0.5,0.5]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0.,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Reflection_x_offset0(self):
        t=Reflection([1.,0.,0.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([-1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([-1.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_x_offset2(self):
        t=Reflection([-2.,0.,0.],offset=-4)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([3.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([4,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([4,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([3.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_x_offset2_vector(self):
        t=Reflection([1.,0.,0.],offset=[2,0,0])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([3.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([4,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([4,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([3.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_y_offset0(self):
        t=Reflection([0.,1.,0.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,-1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([1.,-2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_y_offset2(self):
        t=Reflection([0.,-2.,0.],offset=-4)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,4,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,3,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,4,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([1.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_y_offset2_vector(self):
        t=Reflection([0.,1.,0.],offset=[0,2,0])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,4,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,3,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,4,1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([1.,2,3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_z_offset0(self):
        t=Reflection([0.,0.,1.])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,1,0]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,-1]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([1.,2,-3]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_z_offset2(self):
        t=Reflection([0.,0.,-2.],offset=-4)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,4.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,1,4]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,3]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([1.,2,1]))<self.ABS_TOL,"s is wrong.")
   def test_Reflection_z_offset2_vector(self):
        t=Reflection([0.,0.,1.],offset=[0,0,2])
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,4.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0,1,4]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0,0,3]))<self.ABS_TOL,"s2 is wrong.")
        s=t([1,2,3])
        self.assertTrue(isinstance(s,numpy.ndarray),"s is not an ndarray object.")
        self.assertTrue(self.__distance(s,numpy.array([1.,2,1]))<self.ABS_TOL,"s is wrong.")
   def test_Rotatation_x_90_0(self):
        t=Rotatation(axis=[1.,0.,0.],point=[1.,0.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,0,1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,-1.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_x_30_0(self):
        t=Rotatation(axis=[1.,0.,0.],point=[1.,0.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,0]),numpy.array([1.,0.,0.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([1.,0.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_x_330_0(self):
        t=Rotatation(axis=[1.,0.,0.],point=[1.,0.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]-math.cos(330*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,0]),numpy.array([1.,0.,0.]))>0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([1.,0.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_x_90(self):
        t=Rotatation(axis=[-1.,0.,0.],point=[2.,0.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,0,-1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,1.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_x_30(self):
        t=Rotatation(axis=[-1.,0.,0.],point=[1.,0.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,0]),numpy.array([-1.,0.,0.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([-1.,0.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_x_330(self):
        t=Rotatation(axis=[-1.,0.,0.],point=[1.,0.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]-math.cos(330*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,0]),numpy.array([-1.,0.,0.]))>0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([-1.,0.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_y_90_0(self):
        t=Rotatation(axis=[0.,1.,0.],point=[0.,1.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([0.,0,-1.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,5,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([1,0.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_y_30_0(self):
        t=Rotatation(axis=[0.,1.,0.],point=[0.,1.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,1.,0.]))<0.,"s0 has wrong orientation.")
        s1=t([0,5,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,5,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([0.,1.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_y_330_0(self):
        t=Rotatation(axis=[0.,1.,0.],point=[0.,1.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,1.,0.]))>0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([0.,1.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_y_90(self):
        t=Rotatation(axis=[0.,-1.,0.],point=[0.,2.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([0.,0,1.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,5,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([-1,0.,0.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_y_30(self):
        t=Rotatation(axis=[0.,-1.,0.],point=[0.,2.,0.],angle=30*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,-1.,0.]))<0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(30*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([0.,-1.,0.]))<0.,"s2 has wrong orientation.")
   def test_Rotatation_y_330(self):
        t=Rotatation(axis=[0.,-1.,0.],point=[0.,2.,0.],angle=330*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,-1.,0.]))>0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-1.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(s2[2]-math.cos(330*DEG))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[0,0,1]),numpy.array([0.,-1.,0.]))>0.,"s2 has wrong orientation.")
   def test_Rotatation_z_90_0(self):
        t=Rotatation(axis=[0.,0.,1.],point=[0.,0.,1.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([0.,1,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([-5.,0,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_30_0(self):
        t=Rotatation(axis=[0.,0.,1.],point=[0.,0.,1.],angle=30*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,0.,1.]))<0.,"s0 has wrong orientation.")
        s1=t([0,5,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-5.**2)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]/5.-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,5,0]),numpy.array([0.,0.,1.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_330_0(self):
        t=Rotatation(axis=[0.,0.,1.],point=[0.,0.,1.],angle=330*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,0.,1.]))>0.,"s0 has wrong orientation.")
        s1=t([0,5,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-5.**2)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]/5.-math.cos(330*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,0]),numpy.array([0.,0.,1.]))>0.,"s1 has wrong orientation.")
   def test_Rotatation_z_90(self):
        t=Rotatation(axis=[0.,0.,-1.],point=[0.,0.,2.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([0.,-1,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,5,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([5.,0,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_30(self):
        t=Rotatation(axis=[0.,0.,-1.],point=[0.,0.,2.],angle=30*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(30*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,0.,-1.]))<0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,0]),numpy.array([0.,0.,-1.]))<0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_330(self):
        t=Rotatation(axis=[0.,0.,-1.],point=[0.,0.,2.],angle=330*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-1.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(s0[0]-math.cos(330*DEG))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,0,0]),numpy.array([0.,0.,-1.]))>0.,"s0 has wrong orientation.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-1.)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(s1[1]-math.cos(30*DEG))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,0]),numpy.array([0.,0.,-1.]))>0.,"s1 has wrong orientation.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_x_90_1(self):
        t=Rotatation(point=[0.,0.,1.],axis=[1.,0.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,1,1.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,1,2.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([0.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_y_90_1(self):
        t=Rotatation(point=[1.,0.,0.],axis=[0.,1.,0.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,0,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([1.,1,1.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([2.,0,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_z_90_1(self):
        t=Rotatation(point=[0.,1.,0.],axis=[0.,0.,1.],angle=90*DEG)
        s0=t([1,0,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(self.__distance(s0,numpy.array([1.,2,0.]))<self.ABS_TOL,"s0 is wrong.")
        s1=t([0,1,0])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(self.__distance(s1,numpy.array([0.,1,0.]))<self.ABS_TOL,"s1 is wrong.")
        s2=t([0,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(self.__distance(s2,numpy.array([1.,1,1.]))<self.ABS_TOL,"s2 is wrong.")
   def test_Rotatation_diag_90_0(self):
        t=Rotatation(axis=[1.,1.,1.],angle=90*DEG)
        s0=t([1,-1,0])
        self.assertTrue(isinstance(s0,numpy.ndarray),"s0 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s0,s0)-2.)<self.ABS_TOL,"s0 length is wrong.")
        self.assertTrue(abs(numpy.dot(s0,numpy.array([1,-1,0])))<self.ABS_TOL,"s0 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s0,[1,-1,0]),numpy.array([1.,1.,1.]))<0.,"s0 has wrong orientation.")
        s1=t([0,1,-1])
        self.assertTrue(isinstance(s1,numpy.ndarray),"s1 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s1,s1)-2.)<self.ABS_TOL,"s1 length is wrong.")
        self.assertTrue(abs(numpy.dot(s1,numpy.array([0,1,-1])))<self.ABS_TOL,"s1 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s1,[0,1,-1]),numpy.array([1.,1.,1.]))<0.,"s1 has wrong orientation.")
        s2=t([-1,0,1])
        self.assertTrue(isinstance(s2,numpy.ndarray),"s2 is not an ndarray object.")
        self.assertTrue(abs(numpy.dot(s2,s2)-2.)<self.ABS_TOL,"s2 length is wrong.")
        self.assertTrue(abs(numpy.dot(s2,numpy.array([-1,0,1])))<self.ABS_TOL,"s2 angle is wrong.")
        self.assertTrue(numpy.dot(_cross(s2,[-1,0,1]),numpy.array([1.,1.,1.]))<0.,"s2 has wrong orientation.")
        s3=t([1,1,1])
        self.assertTrue(isinstance(s3,numpy.ndarray),"s3 is not an ndarray object.")
        self.assertTrue(self.__distance(s3,numpy.array([1.,1,1.]))<self.ABS_TOL,"s3 is wrong.")

class Test_PyCAD_Primitives(unittest.TestCase):
   def setUp(self):
         resetGlobalPrimitiveIdCounter()

   def test_Primitive(self):
         p=Primitive()

         id=p.getID()
         self.assertTrue(isinstance(id,int),"id number is not an integer")
         self.assertTrue(not id==Primitive().getID(),"id number is not unique")

         self.assertTrue(p==p.getUnderlyingPrimitive(),"getUnderlyingPrimitive does not return self.")

   def test_ReversePrimitive(self):
         p=Primitive()
         rp=ReversePrimitive(p)
         self.assertTrue(p.getID()==rp.getID(),"reverse primitive does not have same id like source")
         self.assertTrue(p==rp.getUnderlyingPrimitive(),"getUnderlyingPrimitive does return source.")
         self.assertTrue(p == -rp,"reverse or reverse does not return source.")

   def test_Point(self):
       p=Point(1.,2.,3.,local_scale=9.)

       id=p.getID()
       self.assertTrue(isinstance(id,int),"id number is not an integer")
       self.assertTrue(not id==Primitive().getID(),"id number is not unique")

       # check reverse point
       self.assertTrue(p == -p,"reverse is not working.")

       # check history:
       hs=p.getPrimitives()
       self.assertTrue(len(hs)==1,"history must have length 1.")
       self.assertTrue(p in hs,"history must contain point p")

       # check incolved points:
       ps=p.getConstructionPoints()
       self.assertTrue(len(ps)==1,"point set must have length 1.")
       self.assertTrue(p in ps,"point set must contain point p")

       # check coordinates:
       c=p.getCoordinates()
       self.assertTrue(isinstance(c,numpy.ndarray),"coordinates are not an ndarray object.")
       self.assertTrue(c[0]==1.,"x coordinate is not 1.")
       self.assertTrue(c[1]==2.,"y coordinate is not 2.")
       self.assertTrue(c[2]==3.,"z coordinate is not 3.")

       # reset coordinates:
       p.setCoordinates([-1.,-2.,-3.])
       c=p.getCoordinates()
       self.assertTrue(isinstance(c,numpy.ndarray),"new coordinates are not an ndarray object.")
       self.assertTrue(c[0]==-1.,"new x coordinate is not -1.")
       self.assertTrue(c[1]==-2.,"new y coordinate is not -2.")
       self.assertTrue(c[2]==-3.,"new z coordinate is not -3.")

       # check for a colocated point:
       self.assertTrue(p.isColocated(Point(-1.,-2.,-3.)),"colocation not detected.")
       self.assertTrue(not p.isColocated(numpy.array([-1.,-2.,-3.])),"colocation with ndarray representation not detected.")
       self.assertTrue(not p.isColocated(Point(1.,-2.,-3.)),"false colocation detected.")
       self.assertTrue(not p.isColocated(Point(0.,0.,0.)),"false colocation with origin detected.")

       # check for local length scale
       l=p.getLocalScale()
       self.assertTrue(l==9.,"refinement scale is not 9.")

       # check for new local length scale
       p.setLocalScale(3.)
       l=p.getLocalScale()
       self.assertTrue(l==3.,"new refinement scale is not 3.")

       # negative value shouldn't work.
       self.assertRaises(ValueError,p.setLocalScale,-3.)

       # copy:
       an_other_p=p.copy()
       self.assertTrue(isinstance(an_other_p ,Point),"copy is not a point")
       self.assertTrue(not an_other_p.getID() == p.getID(),"copy has same Id")
       self.assertTrue(p.isColocated(an_other_p),"p is not colocated with its copy.")
       self.assertTrue(an_other_p.isColocated(p),"the copy is not colocated with p.")
       self.assertTrue(an_other_p.getLocalScale()==3.,"copy has wrong local scale.")

       # modify by Transformation:
       p.modifyBy(Dilation(-1))
       self.assertTrue(p.isColocated(Point(1.,2.,3.)),"in-place transformation failed")

       # apply Transformation:
       dil_p=p.apply(Dilation(4))
       self.assertTrue(dil_p.isColocated(Point(4.,8.,12.)),"applying transformation failed")
       self.assertTrue(not dil_p.getID() == p.getID(),"transformed point has same Id")
       self.assertTrue(dil_p.getLocalScale()==3.,"transformed point  has wrong local scale.")

       # overloaded add:
       shift_p=p+[1,1,1]
       self.assertTrue(shift_p.isColocated(Point(2,3.,4)),"applying shift by list failed")
       self.assertTrue(not shift_p.getID() == p.getID(),"shift by list has same Id")
       self.assertTrue(shift_p.getLocalScale()==3.,"shift by list has wrong local scale.")

       shift_p=p+numpy.array([1,1,1])
       self.assertTrue(shift_p.isColocated(Point(2,3.,4)),"applying shift by ndarray failed")
       self.assertTrue(not shift_p.getID() == p.getID(),"shift by ndarray has same Id")
       self.assertTrue(shift_p.getLocalScale()==3.,"shift by ndarray has wrong local scale.")
       # overloaded minus
       shift_p=p-[1,1,1]
       self.assertTrue(shift_p.isColocated(Point(0,1,2.)),"applying shift by -list failed")
       self.assertTrue(not shift_p.getID() == p.getID(),"shift by -list has same Id")
       self.assertTrue(shift_p.getLocalScale()==3.,"shift by -list has wrong local scale.")

       shift_p=p-numpy.array([1,1,1])
       self.assertTrue(shift_p.isColocated(Point(0,1,2.)),"applying shift by -ndarray failed")
       self.assertTrue(not shift_p.getID() == p.getID(),"shift by -ndarray has same Id")
       self.assertTrue(shift_p.getLocalScale()==3.,"shift by -ndarray has wrong local scale.")
       # overloaded inplace add:
       p+=[1,1,1]
       self.assertTrue(p.isColocated(Point(2,3.,4)),"modification by list shift failed")

       p+=numpy.array([1,1,1])
       self.assertTrue(p.isColocated(Point(3,4,5)),"modification by ndarray shift failed")

       # overloaded inplace add:
       p-=[1,1,1]
       self.assertTrue(p.isColocated(Point(2,3,4)),"modification by -list shift failed")

       p-=numpy.array([1,1,1])
       self.assertTrue(p.isColocated(Point(1,2.,3)),"modification by -ndarray shift failed")

       #overloaded multiplication:
       mult_p=2*p
       self.assertTrue(mult_p.isColocated(Point(2,4,6)),"applying int factor failed")
       self.assertTrue(not mult_p.getID() == p.getID(),"shift by int factor has same Id")
       self.assertTrue(mult_p.getLocalScale()==3.,"shift by int factor has wrong local scale.")

       mult_p=2.*p
       self.assertTrue(mult_p.isColocated(Point(2,4,6)),"applying float factor failed")
       self.assertTrue(not mult_p.getID() == p.getID(),"shift by float factor has same Id")
       self.assertTrue(mult_p.getLocalScale()==3.,"shift by float factor has wrong local scale.")

       mult_p=Dilation(2)*p
       self.assertTrue(mult_p.isColocated(Point(2,4,6)),"applying Dilation factor failed")
       self.assertTrue(not mult_p.getID() == p.getID(),"shift by Dilation factor has same Id")
       self.assertTrue(mult_p.getLocalScale()==3.,"shift by Dilation factor has wrong local scale.")

       #overloaded inplace multiplication:
       p*=2
       self.assertTrue(p.isColocated(Point(2,4,6)),"applying in-place int factor failed")

       p*=2.
       self.assertTrue(p.isColocated(Point(4,8,12)),"applying in-place float factor failed")

       p*=Dilation(2)
       self.assertTrue(p.isColocated(Point(8,16,24)),"applying in-place Dilation factor failed")

   def test_Spline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)

        self.assertRaises(ValueError,Spline,p0)
        c=Spline(p0,p1,p2,p3)

        self.assertTrue(len(c) == 4, "wrong spline curve length")
        self.assertTrue(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.assertTrue(c.getEndPoint()==p3, "wrong end point of spline curve")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p1),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(Spline(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.assertTrue(not c.isColocated(Spline(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.assertTrue(c.isColocated(Spline(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Spline(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.assertTrue(co[0]==p0, "1st control point is wrong.")
        self.assertTrue(co[1]==p1, "2nd control point is wrong.")
        self.assertTrue(co[2]==p2, "3rd control point is wrong.")
        self.assertTrue(co[3]==p3, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.assertTrue(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.assertTrue(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.assertTrue(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.assertTrue(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 5, "number of primitives in history is wrong.")
        self.assertTrue(p0 in h, "missing p0 in history.")
        self.assertTrue(p1 in h, "missing p1 in history.")
        self.assertTrue(p2 in h, "missing p2 in history.")
        self.assertTrue(p3 in h, "missing p3 in history.")
        self.assertTrue(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.assertTrue(not cp == c, "copy returns same spline curve.")
        self.assertTrue(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.assertTrue(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.assertTrue(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.assertTrue(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.assertTrue(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.assertTrue(c.isColocated(Spline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.assertTrue(p0 == cp[0],"1st new point after Dilation.")
        self.assertTrue(p1 == cp[1],"2nd new point after Dilation.")
        self.assertTrue(p2 == cp[2],"3rd new point after Dilation.")
        self.assertTrue(p3 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.assertTrue(dc.isColocated(Spline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.assertTrue(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.assertTrue(dccp[0].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.assertTrue(dccp[1].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.assertTrue(dccp[2].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p3 == dccp[3],"4th point of Dilation is identical to source.")
        self.assertTrue(dccp[3].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")

   def test_ReverseSpline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        CC0=Spline(p0,p1,p2,p3)
        c=-CC0

        self.assertTrue(len(c) == 4, "wrong reverse spline curve length")
        self.assertTrue(c.getStartPoint()==p3, "wrong start point of reverse spline curve")
        self.assertTrue(c.getEndPoint()==p0, "wrong end point of reverse spline curve")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p1),"reverse spline is colocated with point.")
        self.assertTrue(not c.isColocated(Spline(p0,p1,p2)),"reverse spline is colocated with spline of different length.")
        self.assertTrue(not c.isColocated(Spline(p0,p1,p4,p3)),"reverse spline is colocated with spline with different point.")
        self.assertTrue(c.isColocated(Spline(p0,p1,p2,p3)),"reverse spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(c.isColocated(Spline(p3,p2,p1,p0)),"reverse spline is not colocated with spline with same points.")
        self.assertTrue(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.assertTrue(co[0]==p3, "1st control point is wrong.")
        self.assertTrue(co[1]==p2, "2nd control point is wrong.")
        self.assertTrue(co[2]==p1, "3rd control point is wrong.")
        self.assertTrue(co[3]==p0, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.assertTrue(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.assertTrue(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.assertTrue(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.assertTrue(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 5, "number of primitives in history is wrong.")
        self.assertTrue(p0 in h, "missing p0 in history.")
        self.assertTrue(p1 in h, "missing p1 in history.")
        self.assertTrue(p2 in h, "missing p2 in history.")
        self.assertTrue(p3 in h, "missing p3 in history.")
        self.assertTrue(CC0 in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.assertTrue(not cp == c, "copy returns same spline curve.")
        self.assertTrue(not cp == CC0, "copy returns same spline curve.")
        self.assertTrue(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.assertTrue(not p3 == cpcp[0],"1st point of deep copy and souce are the same.")
        self.assertTrue(not p2 == cpcp[1],"2st point of deep copy and source are the same.")
        self.assertTrue(not p1 == cpcp[2],"3st point of deep copy and source are the same.")
        self.assertTrue(not p0 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.assertTrue(c.isColocated(Spline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.assertTrue(p3 == cp[0],"1st new point after Dilation.")
        self.assertTrue(p2 == cp[1],"2nd new point after Dilation.")
        self.assertTrue(p1 == cp[2],"3rd new point after Dilation.")
        self.assertTrue(p0 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.assertTrue(dc.isColocated(Spline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.assertTrue(dccp[0].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(dccp[1].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(dccp[2].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(dccp[3].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")

   def test_BezierCurve(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        self.assertRaises(ValueError,BezierCurve,p0)
        c=BezierCurve(p0,p1,p2,p3)

        self.assertTrue(len(c) == 4, "wrong spline curve length")
        self.assertTrue(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.assertTrue(c.getEndPoint()==p3, "wrong end point of spline curve")

        self.assertTrue(not c.isColocated(p1),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(BezierCurve(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.assertTrue(not c.isColocated(BezierCurve(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.assertTrue(c.isColocated(BezierCurve(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(BezierCurve(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.assertTrue(co[0]==p0, "1st control point is wrong.")
        self.assertTrue(co[1]==p1, "2nd control point is wrong.")
        self.assertTrue(co[2]==p2, "3rd control point is wrong.")
        self.assertTrue(co[3]==p3, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.assertTrue(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.assertTrue(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.assertTrue(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.assertTrue(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 5, "number of primitives in history is wrong.")
        self.assertTrue(p0 in h, "missing p0 in history.")
        self.assertTrue(p1 in h, "missing p1 in history.")
        self.assertTrue(p2 in h, "missing p2 in history.")
        self.assertTrue(p3 in h, "missing p3 in history.")
        self.assertTrue(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.assertTrue(not cp == c, "copy returns same spline curve.")
        self.assertTrue(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.assertTrue(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.assertTrue(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.assertTrue(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.assertTrue(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.assertTrue(c.isColocated(BezierCurve(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.assertTrue(p0 == cp[0],"1st new point after Dilation.")
        self.assertTrue(p1 == cp[1],"2nd new point after Dilation.")
        self.assertTrue(p2 == cp[2],"3rd new point after Dilation.")
        self.assertTrue(p3 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.assertTrue(dc.isColocated(BezierCurve(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.assertTrue(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.assertTrue(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.assertTrue(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.assertTrue(not p3 == dccp[3],"4th point of Dilation is identical to source.")

   def test_BSpline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        self.assertRaises(ValueError,BSpline,p0)
        c=BSpline(p0,p1,p2,p3)

        self.assertTrue(len(c) == 4, "wrong spline curve length")
        self.assertTrue(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.assertTrue(c.getEndPoint()==p3, "wrong end point of spline curve")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p1),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(BSpline(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.assertTrue(not c.isColocated(BSpline(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.assertTrue(c.isColocated(BSpline(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(BSpline(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.assertTrue(co[0]==p0, "1st control point is wrong.")
        self.assertTrue(co[1]==p1, "2nd control point is wrong.")
        self.assertTrue(co[2]==p2, "3rd control point is wrong.")
        self.assertTrue(co[3]==p3, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.assertTrue(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.assertTrue(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.assertTrue(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.assertTrue(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 5, "number of primitives in history is wrong.")
        self.assertTrue(p0 in h, "missing p0 in history.")
        self.assertTrue(p1 in h, "missing p1 in history.")
        self.assertTrue(p2 in h, "missing p2 in history.")
        self.assertTrue(p3 in h, "missing p3 in history.")
        self.assertTrue(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.assertTrue(not cp == c, "copy returns same spline curve.")
        self.assertTrue(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.assertTrue(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.assertTrue(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.assertTrue(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.assertTrue(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.assertTrue(c.isColocated(BSpline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.assertTrue(p0 == cp[0],"1st new point after Dilation.")
        self.assertTrue(p1 == cp[1],"2nd new point after Dilation.")
        self.assertTrue(p2 == cp[2],"3rd new point after Dilation.")
        self.assertTrue(p3 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.assertTrue(dc.isColocated(BSpline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.assertTrue(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.assertTrue(dccp[0].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.assertTrue(dccp[1].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.assertTrue(dccp[2].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p3 == dccp[3],"4th point of Dilation is identical to source.")
        self.assertTrue(dccp[3].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")

   def test_ReverseBSpline(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p2=Point(2,2,2,0.3)
        p3=Point(3,3,3,0.4)
        p4=Point(1,2,3)
        CC0=BSpline(p0,p1,p2,p3)
        c=-CC0

        self.assertTrue(len(c) == 4, "wrong spline curve length")
        self.assertTrue(c.getStartPoint()==p3, "wrong start point of spline curve")
        self.assertTrue(c.getEndPoint()==p0, "wrong end point of spline curve")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p1),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(BSpline(p0,p1,p2)),"spline is colocated with spline of different length.")
        self.assertTrue(not c.isColocated(BSpline(p0,p1,p4,p3)),"spline is colocated with spline with different point.")
        self.assertTrue(c.isColocated(BSpline(p0,p1,p2,p3)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(BSpline(p3,p2,p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(p0,p1,p2,p3)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.assertTrue(co[0]==p3, "1st control point is wrong.")
        self.assertTrue(co[1]==p2, "2nd control point is wrong.")
        self.assertTrue(co[2]==p1, "3rd control point is wrong.")
        self.assertTrue(co[3]==p0, "4th control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.assertTrue(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.assertTrue(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")
        self.assertTrue(co[2].getLocalScale() == 3., "new local scale of 3rd control point is wrong.")
        self.assertTrue(co[3].getLocalScale() == 3., "new local scale of 4th control point is wrong.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 5, "number of primitives in history is wrong.")
        self.assertTrue(p0 in h, "missing p0 in history.")
        self.assertTrue(p1 in h, "missing p1 in history.")
        self.assertTrue(p2 in h, "missing p2 in history.")
        self.assertTrue(p3 in h, "missing p3 in history.")
        self.assertTrue(CC0 in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.assertTrue(not cp == c, "copy returns same spline curve.")
        self.assertTrue(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.assertTrue(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.assertTrue(not p1 == cpcp[1],"2st point of deep copy and source are the same.")
        self.assertTrue(not p2 == cpcp[2],"3st point of deep copy and source are the same.")
        self.assertTrue(not p3 == cpcp[3],"4st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.assertTrue(c.isColocated(BSpline(Point(0,0,0),Point(-1,-1,-1),Point(-2,-2,-2),Point(-3,-3,-3))),"inplace dilation is wrong.")
        self.assertTrue(p3 == cp[0],"1st new point after Dilation.")
        self.assertTrue(p2 == cp[1],"2nd new point after Dilation.")
        self.assertTrue(p1 == cp[2],"3rd new point after Dilation.")
        self.assertTrue(p0 == cp[3],"4th new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.assertTrue(dc.isColocated(BSpline(Point(0,0,0),Point(1,1,1),Point(2,2,2),Point(3,3,3))),"dilation is wrong.")
        self.assertTrue(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.assertTrue(dccp[0].isColocated(Point(3,3,3)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.assertTrue(dccp[1].isColocated(Point(2,2,2)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p2 == dccp[2],"3rd point of Dilation is identical to source.")
        self.assertTrue(dccp[2].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p3 == dccp[3],"4th point of Dilation is identical to source.")
        self.assertTrue(dccp[3].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")

   def test_LineSegment(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p4=Point(1,2,3)
        self.assertRaises(TypeError,Line,p0)
        self.assertRaises(TypeError,Line,p0,p1,p4)

        c=Line(p0,p1)

        self.assertTrue(len(c) == 2, "wrong spline curve length")
        self.assertTrue(c.getStartPoint()==p0, "wrong start point of spline curve")
        self.assertTrue(c.getEndPoint()==p1, "wrong end point of spline curve")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p1),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(Line(p0,p4)),"spline is colocated with spline with different point.")
        self.assertTrue(c.isColocated(Line(p0,p1)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Line(p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(p0,p1,p4)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.assertTrue(co[0]==p0, "1st control point is wrong.")
        self.assertTrue(co[1]==p1, "2nd control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.assertTrue(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.assertTrue(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 3, "number of primitives in history is wrong.")
        self.assertTrue(p0 in h, "missing p0 in history.")
        self.assertTrue(p1 in h, "missing p1 in history.")
        self.assertTrue(c in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.assertTrue(not cp == c, "copy returns same spline curve.")
        self.assertTrue(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.assertTrue(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.assertTrue(not p1 == cpcp[1],"2st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.assertTrue(c.isColocated(Line(Point(0,0,0),Point(-1,-1,-1))),"inplace dilation is wrong.")
        self.assertTrue(p0 == cp[0],"1st new point after Dilation.")
        self.assertTrue(p1 == cp[1],"2nd new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.assertTrue(dc.isColocated(Line(Point(0,0,0),Point(1,1,1))),"dilation is wrong.")
        self.assertTrue(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.assertTrue(dccp[0].isColocated(Point(0,0,0)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.assertTrue(dccp[1].isColocated(Point(1,1,1)),"2st point of Dilation is is wrongly located.")

        self.assertTrue(dc.getElementDistribution() == None, "element distribution set.")
        dc.setElementDistribution(10,0.2,False)
        d=dc.getElementDistribution()
        self.assertTrue(d[0] == 10, "number of element is wrong.")
        self.assertTrue(d[1] == 0.2, "propagation factor is wrong.")
        self.assertTrue(d[2] == False, "bump flag wrong")
        dc.resetElementDistribution()
        self.assertTrue(dc.getElementDistribution() == None, "resetted element distribution set.")

   def test_ReverseLineSegment(self):
        p0=Point(0,0,0,0.1)
        p1=Point(1,1,1,0.2)
        p4=Point(1,2,3)
        self.assertRaises(TypeError,Line,p0)
        self.assertRaises(TypeError,Line,p0,p1,p4)

        CC0=Line(p0,p1)
        c=-CC0

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(len(c) == 2, "wrong spline curve length")
        self.assertTrue(c.getStartPoint()==p1, "wrong start point of spline curve")
        self.assertTrue(c.getEndPoint()==p0, "wrong end point of spline curve")

        self.assertTrue(not c.isColocated(p1),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(Line(p0,p4)),"spline is colocated with spline with different point.")
        self.assertTrue(c.isColocated(Line(p0,p1)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Line(p1,p0)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(p0,p1,p4)),"spline curve is identified with curve.")

        co=c.getControlPoints()
        self.assertTrue(co[0]==p1, "1st control point is wrong.")
        self.assertTrue(co[1]==p0, "2nd control point is wrong.")

        c.setLocalScale(3.)
        co=c.getControlPoints()
        self.assertTrue(co[0].getLocalScale() == 3., "new local scale of 1st control point is wrong.")
        self.assertTrue(co[1].getLocalScale() == 3., "new local scale of 2nd control point is wrong.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 3, "number of primitives in history is wrong.")
        self.assertTrue(p0 in h, "missing p0 in history.")
        self.assertTrue(p1 in h, "missing p1 in history.")
        self.assertTrue(CC0 in h, "missing spline curve in history.")

        cp=c.copy()
        cpcp=cp.getControlPoints()
        self.assertTrue(not cp == c, "copy returns same spline curve.")
        self.assertTrue(c.isColocated(cp),"spline curve is not colocated with its copy.")
        self.assertTrue(not p0 == cpcp[0],"1st point of deep copy and source are the same.")
        self.assertTrue(not p1 == cpcp[1],"2st point of deep copy and source are the same.")

        c.modifyBy(Dilation(-1.))
        cp=c.getControlPoints()
        self.assertTrue(c.isColocated(Line(Point(0,0,0),Point(-1,-1,-1))),"inplace dilation is wrong.")
        self.assertTrue(p1 == cp[0],"1st new point after Dilation.")
        self.assertTrue(p0 == cp[1],"2nd new point after Dilation.")

        dc=c.apply(Dilation(-1.))
        dccp=dc.getControlPoints()
        self.assertTrue(dc.isColocated(Line(Point(0,0,0),Point(1,1,1))),"dilation is wrong.")
        self.assertTrue(not p0 == dccp[0],"1st point of Dilation is identical to source.")
        self.assertTrue(dccp[0].isColocated(Point(1,1,1)),"1st point of Dilation is is wrongly located.")
        self.assertTrue(not p1 == dccp[1],"2nd point of Dilation is identical to source.")
        self.assertTrue(dccp[1].isColocated(Point(0,0,0)),"2st point of Dilation is is wrongly located.")

        self.assertTrue(dc.getElementDistribution() == None, "element distribution set.")
        dc.setElementDistribution(10,0.2,False)
        d=dc.getElementDistribution()
        self.assertTrue(d[0] == 10, "number of element is wrong.")
        self.assertTrue(d[1] == 0.2, "propagation factor is wrong.")
        self.assertTrue(d[2] == False, "bump flag wrong")
        dc.resetElementDistribution()
        self.assertTrue(dc.getElementDistribution() == None, "resetted element distribution set.")

   def test_Arc(self):
        center=Point(0,0,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
        p4=Point(10,2,3)

        self.assertRaises(TypeError,Arc,Primitive())

        c=Arc(center,p_start,p_end)

        self.assertTrue(c.getCenterPoint()==center, "wrong center point")
        self.assertTrue(c.getStartPoint()==p_start, "wrong start point")
        self.assertTrue(c.getEndPoint()==p_end, "wrong end point")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p4),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(Arc(p4,p_start,p_end)),"spline is colocated with spline with differnt center point.")
        self.assertTrue(not c.isColocated(Arc(center,p4,p_end)),"spline is colocated with spline with differnt start point.")
        self.assertTrue(not c.isColocated(Arc(center,p_start,p4)),"spline is colocated with spline with differnt end point.")
        self.assertTrue(c.isColocated(Arc(center,p_start,p_end)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Arc(center,p_end,p_start)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(center,p_start,p_end)),"spline curve is identified with curve.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 4, "number of primitives in history is wrong.")
        self.assertTrue(center in h, "missing center in history.")
        self.assertTrue(p_start in h, "missing p_start in history.")
        self.assertTrue(p_end in h, "missing p_end in history.")
        self.assertTrue(c in h, "missing spline curve in history.")


        c.setLocalScale(3.)
        self.assertTrue(c.getCenterPoint().getLocalScale() == 3., "new local scale of center point is wrong.")
        self.assertTrue(c.getStartPoint().getLocalScale() == 3., "new local scale of start point is wrong.")
        self.assertTrue(c.getEndPoint().getLocalScale() == 3., "new local scale of end point is wrong.")

        cp=c.copy()
        self.assertTrue(isinstance(cp,Arc), "copy returns is not an arc.")
        self.assertTrue(not cp == c, "copy returns same arc.")
        self.assertTrue(cp.isColocated(Arc(center,p_start,p_end)),"arc is not colocated with its copy.")
        self.assertTrue(not cp.getCenterPoint()==center, "deep copy has same center point like source")
        self.assertTrue(not cp.getStartPoint()==p_start, "deep copy has same start point like source")
        self.assertTrue(not cp.getEndPoint()==p_end, "deep copy has same end point like source")

        c.modifyBy(Dilation(-1.))
        self.assertTrue(c.isColocated(Arc(Point(0,0,0),Point(-1,-1,-1),Point(-1,-2,-3))),"inplace dilation is wrong.")
        self.assertTrue(c.getCenterPoint() == center,"wrong center point after dilation.")
        self.assertTrue(c.getStartPoint() == p_start,"wrong start point after dilation.")
        self.assertTrue(c.getEndPoint() == p_end,"wrong end point after dilation.")

        dc=c.apply(Dilation(-1.))
        self.assertTrue(dc.isColocated(Arc(Point(0,0,0),Point(1,1,1),Point(1,2,3))),"dilation is wrong.")
        self.assertTrue(not dc.getCenterPoint() == center,"center point of dilation is identical to source.")
        self.assertTrue(dc.getCenterPoint().isColocated(Point(0,0,0)),"center point of dilation is wrong.")
        self.assertTrue(not dc.getStartPoint() == p_start,"start point of dilation is identical to source.")
        self.assertTrue(dc.getStartPoint().isColocated(Point(1,1,1)),"start point of dilation is wrong.")
        self.assertTrue(not dc.getEndPoint() == p_end,"end point of dilation is identical to source.")
        self.assertTrue(dc.getEndPoint().isColocated(Point(1,2,3)),"end point of dilation is wrong.")

        self.assertTrue(dc.getElementDistribution() == None, "element distribution set.")
        dc.setElementDistribution(10,0.2,False)
        d=dc.getElementDistribution()
        self.assertTrue(d[0] == 10, "number of element is wrong.")
        self.assertTrue(d[1] == 0.2, "propagation factor is wrong.")
        self.assertTrue(d[2] == False, "bump flag wrong")
        dc.resetElementDistribution()
        self.assertTrue(dc.getElementDistribution() == None, "resetted element distribution set.")

   def test_ReverseArc(self):
        center=Point(0,0,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
        p4=Point(10,2,3)
        self.assertRaises(TypeError,Arc,Primitive())

        CC0=Arc(center,p_start,p_end)
        c=-CC0

        self.assertTrue(c.getCenterPoint()==center, "wrong center point")
        self.assertTrue(c.getStartPoint()==p_end, "wrong start point")
        self.assertTrue(c.getEndPoint()==p_start, "wrong end point")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p4),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(Arc(p4,p_start,p_end)),"spline is colocated with spline with differnt center point.")
        self.assertTrue(not c.isColocated(Arc(center,p4,p_end)),"spline is colocated with spline with differnt start point.")
        self.assertTrue(not c.isColocated(Arc(center,p_start,p4)),"spline is colocated with spline with differnt end point.")
        self.assertTrue(c.isColocated(Arc(center,p_start,p_end)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Arc(center,p_end,p_start)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(center,p_start,p_end)),"spline curve is identified with curve.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 4, "number of primitives in history is wrong.")
        self.assertTrue(center in h, "missing center in history.")
        self.assertTrue(p_start in h, "missing p_start in history.")
        self.assertTrue(p_end in h, "missing p_end in history.")
        self.assertTrue(CC0 in h, "missing spline curve in history.")


        c.setLocalScale(3.)
        self.assertTrue(c.getCenterPoint().getLocalScale() == 3., "new local scale of center point is wrong.")
        self.assertTrue(c.getStartPoint().getLocalScale() == 3., "new local scale of start point is wrong.")
        self.assertTrue(c.getEndPoint().getLocalScale() == 3., "new local scale of end point is wrong.")

        cp=c.copy()
        self.assertTrue(isinstance(cp,ReverseArc), "copy returns is not an arc.")
        self.assertTrue(not cp == c, "copy returns same arc.")
        self.assertTrue(cp.isColocated(Arc(center,p_end,p_start)),"arc is not colocated with its copy.")
        self.assertTrue(not cp.getCenterPoint()==center, "deep copy has same center point like source")
        self.assertTrue(not cp.getStartPoint()==p_start, "deep copy has same start point like source")
        self.assertTrue(not cp.getEndPoint()==p_end, "deep copy has same end point like source")

        c.modifyBy(Dilation(-1.))
        self.assertTrue(c.isColocated(Arc(Point(0,0,0),Point(-1,-1,-1),Point(-1,-2,-3))),"inplace dilation is wrong.")
        self.assertTrue(c.getCenterPoint() == center,"wrong center point after dilation.")
        self.assertTrue(c.getStartPoint() == p_end,"wrong start point after dilation.")
        self.assertTrue(c.getEndPoint() == p_start,"wrong end point after dilation.")

        dc=c.apply(Dilation(-1.))
        self.assertTrue(dc.isColocated(Arc(Point(0,0,0),Point(1,1,1),Point(1,2,3))),"dilation is wrong.")
        self.assertTrue(not dc.getCenterPoint() == center,"center point of dilation is identical to source.")
        self.assertTrue(dc.getCenterPoint().isColocated(Point(0,0,0)),"center point of dilation is wrong.")
        self.assertTrue(not dc.getStartPoint() == p_start,"start point of dilation is identical to source.")
        self.assertTrue(dc.getStartPoint().isColocated(Point(1,2,3)),"start point of dilation is wrong.")
        self.assertTrue(not dc.getEndPoint() == p_end,"end point of dilation is identical to source.")
        self.assertTrue(dc.getEndPoint().isColocated(Point(1,1,1)),"end point of dilation is wrong.")

        self.assertTrue(dc.getElementDistribution() == None, "element distribution set.")
        dc.setElementDistribution(10,0.2,False)
        d=dc.getElementDistribution()
        self.assertTrue(d[0] == 10, "number of element is wrong.")
        self.assertTrue(d[1] == 0.2, "propagation factor is wrong.")
        self.assertTrue(d[2] == False, "bump flag wrong")
        dc.resetElementDistribution()
        self.assertTrue(dc.getElementDistribution() == None, "resetted element distribution set.")

   def test_Ellipse(self):
        center=Point(0,0,0,0.1)
        main_axis_point=Point(0,1,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
        p4=Point(10,2,3)
        self.assertRaises(TypeError,Ellipse,Primitive())
        self.assertRaises(TypeError,Ellipse,center,center,p_start,p_end)
        self.assertRaises(TypeError,Ellipse,center,main_axis_point,p_start,p_start)


        c=Ellipse(center,main_axis_point,p_start,p_end)

        self.assertTrue(c.getCenterPoint()==center, "wrong center point")
        self.assertTrue(c.getStartPoint()==p_start, "wrong start point")
        self.assertTrue(c.getEndPoint()==p_end, "wrong end point")
        self.assertTrue(c.getPointOnMainAxis()==main_axis_point, "wrong point on main axis")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p4),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(Ellipse(center,main_axis_point,p4,p_end)),"spline is colocated with spline with differnt start point.")
        self.assertTrue(not c.isColocated(Ellipse(center,main_axis_point,p_start,p4)),"spline is colocated with spline with differnt end point.")
        self.assertTrue(not c.isColocated(Ellipse(center,p4,p_start,p_end)),"spline is colocated with spline with differnt main axis point.")
        self.assertTrue(not c.isColocated(Ellipse(p4,main_axis_point,p_start,p_end)),"spline is colocated with spline with differnt center.")
        self.assertTrue(c.isColocated(Ellipse(center,main_axis_point,p_start,p_end)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Ellipse(center,main_axis_point,p_start,p_end)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Ellipse(center,main_axis_point,p_end,p_start)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(center,main_axis_point,p_start,p_end)),"spline curve is identified with curve.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 5, "number of primitives in history is wrong.")
        self.assertTrue(center in h, "missing center in history.")
        self.assertTrue(p_start in h, "missing p_start in history.")
        self.assertTrue(p_end in h, "missing p_end in history.")
        self.assertTrue(main_axis_point in h, "missing main_axis_point in history.")
        self.assertTrue(c in h, "missing spline curve in history.")


        c.setLocalScale(3.)
        self.assertTrue(c.getCenterPoint().getLocalScale() == 3., "new local scale of center point is wrong.")
        self.assertTrue(c.getStartPoint().getLocalScale() == 3., "new local scale of start point is wrong.")
        self.assertTrue(c.getEndPoint().getLocalScale() == 3., "new local scale of end point is wrong.")
        self.assertTrue(c.getPointOnMainAxis().getLocalScale() == 3., "new local scale of point on main axis is wrong.")

        cp=c.copy()
        self.assertTrue(isinstance(cp,Ellipse), "copy returns is not an arc.")
        self.assertTrue(not cp == c, "copy returns same arc.")
        self.assertTrue(cp.isColocated(Ellipse(center,main_axis_point,p_start,p_end)),"arc is not colocated with its copy.")
        self.assertTrue(not cp.getCenterPoint()==center, "deep copy has same center point like source")
        self.assertTrue(not cp.getStartPoint()==p_start, "deep copy has same start point like source")
        self.assertTrue(not cp.getEndPoint()==p_end, "deep copy has same end point like source")
        self.assertTrue(not cp.getPointOnMainAxis()==main_axis_point, "deep copy has same point on main axis like source")

        c.modifyBy(Dilation(-1.))
        self.assertTrue(c.isColocated(Ellipse(Point(0,0,0),Point(0,-1,0),Point(-1,-1,-1),Point(-1,-2,-3))),"inplace dilation is wrong.")
        self.assertTrue(c.getCenterPoint() == center,"wrong center point after dilation.")
        self.assertTrue(c.getStartPoint() == p_start,"wrong start point after dilation.")
        self.assertTrue(c.getEndPoint() == p_end,"wrong end point after dilation.")
        self.assertTrue(c.getPointOnMainAxis() == main_axis_point,"wrong point on main axis  after dilation.")

        #=====================
        dc=c.apply(Dilation(-1.))
        self.assertTrue(dc.isColocated(Ellipse(Point(0,0,0),Point(0,1,0),Point(1,1,1),Point(1,2,3))),"dilation is wrong.")
        self.assertTrue(not dc.getCenterPoint() == center,"center point of dilation is identical to source.")
        self.assertTrue(dc.getCenterPoint().isColocated(Point(0,0,0)),"center point of dilation is wrong.")
        self.assertTrue(not dc.getStartPoint() == p_start,"start point of dilation is identical to source.")
        self.assertTrue(dc.getStartPoint().isColocated(Point(1,1,1)),"start point of dilation is wrong.")
        self.assertTrue(not dc.getEndPoint() == p_end,"end point of dilation is identical to source.")
        self.assertTrue(dc.getEndPoint().isColocated(Point(1,2,3)),"end point of dilation is wrong.")
        self.assertTrue(not dc.getPointOnMainAxis() == main_axis_point,"point on main axis is identical to source.")
        self.assertTrue(dc.getPointOnMainAxis().isColocated(Point(0,1,0)),"point on main axis of dilation is wrong.")

        self.assertTrue(dc.getElementDistribution() == None, "element distribution set.")
        dc.setElementDistribution(10,0.2,False)
        d=dc.getElementDistribution()
        self.assertTrue(d[0] == 10, "number of element is wrong.")
        self.assertTrue(d[1] == 0.2, "propagation factor is wrong.")
        self.assertTrue(d[2] == False, "bump flag wrong")
        dc.resetElementDistribution()
        self.assertTrue(dc.getElementDistribution() == None, "resetted element distribution set.")

   def test_ReverseEllipse(self):
        center=Point(0,0,0,0.1)
        main_axis_point=Point(0,1,0,0.1)
        p_start=Point(1,1,1,0.2)
        p_end=Point(1,2,3)
        p4=Point(10,2,3)
        self.assertRaises(TypeError,Ellipse,Primitive())

        CC0=Ellipse(center,main_axis_point,p_start,p_end)
        c=-CC0

        self.assertTrue(c.getCenterPoint()==center, "wrong center point")
        self.assertTrue(c.getStartPoint()==p_end, "wrong start point")
        self.assertTrue(c.getEndPoint()==p_start, "wrong end point")
        self.assertTrue(c.getPointOnMainAxis()==main_axis_point, "wrong point on main axis")

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p4),"spline is colocated with point.")
        self.assertTrue(not c.isColocated(Ellipse(center,main_axis_point,p4,p_start)),"spline is colocated with spline with differnt start point.")
        self.assertTrue(not c.isColocated(Ellipse(center,main_axis_point,p_end,p4)),"spline is colocated with spline with differnt end point.")
        self.assertTrue(not c.isColocated(Ellipse(center,p4,p_end,p_start)),"spline is colocated with spline with differnt main axis point.")
        self.assertTrue(not c.isColocated(Ellipse(p4,main_axis_point,p_end,p_start)),"spline is colocated with spline with differnt center.")
        self.assertTrue(c.isColocated(Ellipse(center,main_axis_point,p_end,p_start)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Ellipse(center,main_axis_point,p_end,p_start)),"spline is not colocated with spline with same points.")
        self.assertTrue(c.isColocated(Ellipse(center,main_axis_point,p_start,p_end)),"spline is not colocated with spline with same points but opposite direction.")
        self.assertTrue(not c.isColocated(Curve(center,main_axis_point,p_start,p_end)),"spline curve is identified with curve.")

        h=c.getPrimitives()
        self.assertTrue(len(h) == 5, "number of primitives in history is wrong.")
        self.assertTrue(center in h, "missing center in history.")
        self.assertTrue(p_start in h, "missing p_start in history.")
        self.assertTrue(p_end in h, "missing p_end in history.")
        self.assertTrue(main_axis_point in h, "missing main_axis_point in history.")
        self.assertTrue(CC0 in h, "missing spline curve in history.")


        c.setLocalScale(3.)
        self.assertTrue(c.getCenterPoint().getLocalScale() == 3., "new local scale of center point is wrong.")
        self.assertTrue(c.getStartPoint().getLocalScale() == 3., "new local scale of start point is wrong.")
        self.assertTrue(c.getEndPoint().getLocalScale() == 3., "new local scale of end point is wrong.")
        self.assertTrue(c.getPointOnMainAxis().getLocalScale() == 3., "new local scale of point on main axis is wrong.")

        cp=c.copy()
        self.assertTrue(isinstance(cp,ReverseEllipse), "copy returns is not an arc.")
        self.assertTrue(not cp == c, "copy returns same arc.")
        self.assertTrue(cp.isColocated(Ellipse(center,main_axis_point,p_end,p_start)),"arc is not colocated with its copy.")
        self.assertTrue(not cp.getCenterPoint()==center, "deep copy has same center point like source")
        self.assertTrue(not cp.getStartPoint()==p_start, "deep copy has same start point like source")
        self.assertTrue(not cp.getEndPoint()==p_end, "deep copy has same end point like source")
        self.assertTrue(not cp.getPointOnMainAxis()==main_axis_point, "deep copy has same point on main axis like source")

        c.modifyBy(Dilation(-1.))
        self.assertTrue(c.isColocated(Ellipse(Point(0,0,0),Point(0,-1,0),Point(-1,-1,-1),Point(-1,-2,-3))),"inplace dilation is wrong.")
        self.assertTrue(c.getCenterPoint() == center,"wrong center point after dilation.")
        self.assertTrue(c.getStartPoint() == p_end,"wrong start point after dilation.")
        self.assertTrue(c.getEndPoint() == p_start,"wrong end point after dilation.")
        self.assertTrue(c.getPointOnMainAxis() == main_axis_point,"wrong point on main axis  after dilation.")

        dc=c.apply(Dilation(-1.))
        self.assertTrue(dc.isColocated(Ellipse(Point(0,0,0),Point(0,1,0),Point(1,1,1),Point(1,2,3))),"dilation is wrong.")
        self.assertTrue(not dc.getCenterPoint() == center,"center point of dilation is identical to source.")
        self.assertTrue(dc.getCenterPoint().isColocated(Point(0,0,0)),"center point of dilation is wrong.")
        self.assertTrue(not dc.getStartPoint() == p_start,"start point of dilation is identical to source.")
        self.assertTrue(dc.getStartPoint().isColocated(Point(1,2,3)),"start point of dilation is wrong.")
        self.assertTrue(not dc.getEndPoint() == p_end,"end point of dilation is identical to source.")
        self.assertTrue(dc.getEndPoint().isColocated(Point(1,1,1)),"end point of dilation is wrong.")
        self.assertTrue(not dc.getPointOnMainAxis() == main_axis_point,"point on main axis is identical to source.")
        self.assertTrue(dc.getPointOnMainAxis().isColocated(Point(0,1,0)),"point on main axis of dilation is wrong.")

        self.assertTrue(dc.getElementDistribution() == None, "element distribution set.")
        dc.setElementDistribution(10,0.2,False)
        d=dc.getElementDistribution()
        self.assertTrue(d[0] == 10, "number of element is wrong.")
        self.assertTrue(d[1] == 0.2, "propagation factor is wrong.")
        self.assertTrue(d[2] == False, "bump flag wrong")
        dc.resetElementDistribution()
        self.assertTrue(dc.getElementDistribution() == None, "resetted element distribution set.")

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
        # self.assertRaises(ValueError,CurveLoop,l01,lx,l20)
        # self.assertRaises(ValueError,CurveLoop,l01,l20,l20)
        # self.assertRaises(ValueError,CurveLoop,l01,l20,ly)

        c=CurveLoop(l01,l20,l12)
        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p4),"CurveLoop is colocated with point.")
        self.assertTrue(c.isColocated(c),"CurveLoop is not colocated with its self.")
        self.assertTrue(c.isColocated(CurveLoop(l01,l12,l20)),"CurveLoop is not colocated with its copy.")
        self.assertTrue(c.isColocated(CurveLoop(l20,l01,l12)),"CurveLoop is not colocated with its copy with shifted points.")
        self.assertTrue(c.isColocated(CurveLoop(l20,l12,l01)),"CurveLoop is not colocated with its copy with shuffled points.")
        self.assertTrue(not c.isColocated(CurveLoop(lx,ly,l12)),"CurveLoop is colocated with different CurveLoop.")

        self.assertTrue(len(c) == 3, "wrong length")

        c.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")


        cc=c.getCurves()
        self.assertTrue(len(cc) == 3, "too many curves.")
        self.assertTrue(l01 in cc, "l01 is missing")
        self.assertTrue(l12 in cc, "l12 is missing")
        self.assertTrue(l20 in cc, "l20 is missing")

        p=c.getPrimitives()
        self.assertTrue(len(p) == 9, "too many primitives.")
        self.assertTrue(l01 in p, "l01 is missing")
        self.assertTrue(l12 in p, "l21 is missing")
        self.assertTrue(l20 in p, "l20 is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")

        cp=c.copy()
        self.assertTrue(isinstance(cp,CurveLoop), "copy returns is not an arc.")
        self.assertTrue(not cp == c, "copy equals source")
        self.assertTrue(cp.isColocated(c),"copy is not colocated with its source.")
        cc=cp.getCurves()
        self.assertTrue(len(cc) == 3, "too many primitives in copy.")
        self.assertTrue(not l01 in cc,"copy uses l01.")
        self.assertTrue(not l12 in cc,"copy uses l12.")
        self.assertTrue(not l20 in cc,"copy uses l20.")

        p0_m=Point(0,0,0)
        p1_m=Point(-1,-1,-1)
        p2_m=Point(-2,-2,-2)
        p3_m=Point(-3,-3,-3)
        p4_m=Point(-1,-2,-3)

        l01_m=Line(p0_m,p1_m)
        l12_m=Arc(p3_m,p1_m,p2_m)
        l20_m=Spline(p2_m,p4_m,p0_m)

        dc=c.apply(Dilation(-1.))
        self.assertTrue(dc.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"dilation is wrong.")
        cc=dc.getCurves()
        self.assertTrue(len(cc) == 3, "too many primitives in dilation result.")
        self.assertTrue(not l01 in cc,"l01 is in dilation result.")
        self.assertTrue(not l12 in cc,"l12 is in dilation result.")
        self.assertTrue(not l20 in cc,"l20 is in dilation result.")

        c.modifyBy(Dilation(-1.))
        self.assertTrue(c.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"inplace dilation is wrong.")
        cc=c.getCurves()
        self.assertTrue(len(cc) == 3, "too many primitives in modified object.")
        self.assertTrue(l01 in cc,"l01 missed in  modified object.")
        self.assertTrue(cc[cc.index(l01)].hasSameOrientation(l01),"l01 in modified object has wrong orientation.")
        self.assertTrue(l12 in cc,"l12 missed in  modified object.")
        self.assertTrue(cc[cc.index(l12)].hasSameOrientation(l12),"l12 in modified object has wrong orientation.")
        self.assertTrue(l20 in cc,"l20 missed in  modified object.")
        self.assertTrue(cc[cc.index(l20)].hasSameOrientation(l20),"l20 in modified object has wrong orientation.")

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

        self.assertTrue(c.hasSameOrientation(c),"has not same orientation like itself")
        self.assertTrue(not c.hasSameOrientation(-c),"has same orientation like -itself")

        self.assertTrue(not c.isColocated(p4),"-CurveLoop is colocated with point.")
        self.assertTrue(c.isColocated(c),"-CurveLoop is not colocated with its self.")
        self.assertTrue(c.isColocated(CurveLoop(l01,l12,l20)),"-CurveLoop is not colocated with its copy.")
        self.assertTrue(c.isColocated(CurveLoop(l20,l01,l12)),"-CurveLoop is not colocated with its copy with shifted points.")
        self.assertTrue(c.isColocated(CurveLoop(l20,l12,l01)),"-CurveLoop is not colocated with its copy with shuffled points.")
        self.assertTrue(not c.isColocated(CurveLoop(lx,ly,l12)),"-CurveLoop is colocated with different CurveLoop.")

        self.assertTrue(len(c) == 3, "wrong length")

        c.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")


        cc=c.getCurves()
        self.assertTrue(len(cc) == 3, "too many curves.")
        self.assertTrue(l01 in cc, "l01 is missing")
        self.assertTrue(l12 in cc, "l12 is missing")
        self.assertTrue(l20 in cc, "l20 is missing")

        p=c.getPrimitives()
        self.assertTrue(len(p) == 9, "too many primitives.")
        self.assertTrue(l01 in p, "l01 is missing")
        self.assertTrue(l12 in p, "l21 is missing")
        self.assertTrue(l20 in p, "l20 is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")

        cp=c.copy()
        self.assertTrue(isinstance(cp,ReverseCurveLoop), "copy returns is not an ReverseCurveLoop.")
        self.assertTrue(not cp == c, "copy equals source")
        self.assertTrue(cp.isColocated(c),"copy is not colocated with its source.")
        cc=cp.getCurves()
        self.assertTrue(len(cc) == 3, "too many primitives in copy.")
        self.assertTrue(not l01 in cc,"copy uses l01.")
        self.assertTrue(not l12 in cc,"copy uses l12.")
        self.assertTrue(not l20 in cc,"copy uses l20.")

        p0_m=Point(0,0,0)
        p1_m=Point(-1,-1,-1)
        p2_m=Point(-2,-2,-2)
        p3_m=Point(-3,-3,-3)
        p4_m=Point(-1,-2,-3)

        l01_m=Line(p0_m,p1_m)
        l12_m=Arc(p3_m,p1_m,p2_m)
        l20_m=Spline(p2_m,p4_m,p0_m)

        dc=c.apply(Dilation(-1.))
        self.assertTrue(dc.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"dilation is wrong.")
        cc=dc.getCurves()
        self.assertTrue(len(cc) == 3, "too many primitives in dilation result.")
        self.assertTrue(not l01 in cc,"l01 is in dilation result.")
        self.assertTrue(not l12 in cc,"l12 is in dilation result.")
        self.assertTrue(not l20 in cc,"l20 is in dilation result.")

        c.modifyBy(Dilation(-1.))
        self.assertTrue(c.isColocated(CurveLoop(l01_m,l12_m,l20_m)),"inplace dilation is wrong.")
        cc=c.getCurves()
        self.assertTrue(len(cc) == 3, "too many primitives in modified object.")
        self.assertTrue(l01 in cc,"l01 missed in  modified object.")
        self.assertTrue(cc[cc.index(l01)].hasSameOrientation(-l01),"l01 in modified object has wrong orientation.")
        self.assertTrue(l12 in cc,"l12 missed in  modified object.")
        self.assertTrue(cc[cc.index(l12)].hasSameOrientation(-l12),"l12 in modified object has wrong orientation.")
        self.assertTrue(l20 in cc,"l20 missed in  modified object.")
        self.assertTrue(cc[cc.index(l20)].hasSameOrientation(-l20),"l20 in modified object has wrong orientation.")

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

        self.assertRaises(TypeError,RuledSurface,l01)

        s=RuledSurface(cl1)

        cl=s.getBoundaryLoop()
        self.assertTrue(cl == cl1, " wrong boundary loops")
        self.assertTrue(cl.hasSameOrientation(cl1),"cl1 has incorrect orientation.")

        self.assertTrue(s.hasSameOrientation(s),"has not same orientation like itself")
        self.assertTrue(not s.hasSameOrientation(-s),"has same orientation like -itself")

        crvs=s.getBoundary()
        self.assertTrue(len(crvs) == 3, "too many boundary corves.")
        self.assertTrue(l01 in crvs, "l01 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l01)].hasSameOrientation(l01),"l01 has incorrect orientation.")
        self.assertTrue(l12_1 in crvs, "l21 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l12_1)].hasSameOrientation(l12_1),"l12_1 has incorrect orientation.")
        self.assertTrue(l20 in crvs, "l20 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l20)].hasSameOrientation(l20),"l12_1 has incorrect orientation.")

        self.assertTrue(not s.isColocated(p4),"RuledSurface is colocated with point.")
        self.assertTrue(s.isColocated(s),"RuledSurface is not colocated with its self.")
        self.assertTrue(s.isColocated(RuledSurface(cl1)),"RuledSurface is not colocated with its copy.")
        self.assertTrue(not s.isColocated(RuledSurface(cl2)),"RuledSurface is colocated with different length")
        self.assertTrue(not s.isColocated(RuledSurface(cl3)),"RuledSurface is colocated with same length.")

        s.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")

        p=s.getPrimitives()
        self.assertTrue(len(p) == 10, "too many primitives.")
        self.assertTrue(cl1 in p, "cl1 is missing")
        self.assertTrue(l01 in p, "l01 is missing")
        self.assertTrue(l12_1 in p, "l21 is missing")
        self.assertTrue(l20 in p, "l20 is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")

        sp=s.copy()
        self.assertTrue(isinstance(sp,RuledSurface), "copy returns is not a RuledSurface.")
        self.assertTrue(not sp == s, "copy equals source")
        self.assertTrue(sp.isColocated(s),"copy is not colocated with its source.")
        cbl=sp.getBoundaryLoop()
        self.assertTrue(not cbl == cl1,"copy uses cl1.")
        cp=sp.getPrimitives()
        self.assertTrue(len(cp) == 10, "copy as too many primitives.")
        self.assertTrue(not cl1 in cp, "copy is using cl1")
        self.assertTrue(not l01 in cp, "copy is using l01")
        self.assertTrue(not l12_1 in cp, "copy is using l21")
        self.assertTrue(not l20 in cp, "copy is using l20")
        self.assertTrue(not p0 in cp, "copy is using p0")
        self.assertTrue(not p1 in cp, "copy is using p1")
        self.assertTrue(not p2 in cp, "copy is using p2")
        self.assertTrue(not p3 in cp, "copy is using p3")
        self.assertTrue(not p4 in cp, "copy is using p4")
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
        self.assertTrue(ds.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"dilation is wrong.")
        cbl=ds.getBoundaryLoop()
        self.assertTrue(not cbl == cl1,"dilation uses cl1.")
        cp=ds.getPrimitives()
        self.assertTrue(len(cp) == 10, "dilation as too many primitives.")
        self.assertTrue(not cl1 in cp, "dilation is using cl1")
        self.assertTrue(not l01 in cp, "dilation is using l01")
        self.assertTrue(not l12_1 in cp, "dilation is using l21")
        self.assertTrue(not l20 in cp, "dilation is using l20")
        self.assertTrue(not p0 in cp, "dilation is using p0")
        self.assertTrue(not p1 in cp, "dilation is using p1")
        self.assertTrue(not p2 in cp, "dilation is using p2")
        self.assertTrue(not p3 in cp, "dilation is using p3")
        self.assertTrue(not p4 in cp, "dilation is using p4")

        s.modifyBy(Dilation(-1.))
        self.assertTrue(s.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"inplace dilation is wrong.")

        p=s.getPrimitives()
        self.assertTrue(len(p) == 10, "inplace dilation has too many primitives.")
        self.assertTrue(cl1 in p, "inplace dilation cl1 is missing")
        self.assertTrue(l01 in p, "inplace dilation l01 is missing")
        self.assertTrue(l12_1 in p, "inplace dilation l21 is missing")
        self.assertTrue(l20 in p, "inplace dilation l20 is missing")
        self.assertTrue(p0 in p, "inplace dilation p0 is missing")
        self.assertTrue(p1 in p, "inplace dilation p1 is missing")
        self.assertTrue(p2 in p, "inplace dilation p2 is missing")
        self.assertTrue(p3 in p, "inplace dilation p3 is missing")
        self.assertTrue(p4 in p, "inplace dilation p4 is missing")

        p=s.getBoundary()
        self.assertTrue(len(p) == 3, "inplace dilation has too many boundary curves.")
        self.assertTrue(l01 in p, "inplace dilation l01 is missing in boundary curves.")
        self.assertTrue(p[p.index(l01)].hasSameOrientation(l01),"l01 in getBoundary after dilation has incorrect orientation.")
        self.assertTrue(l12_1 in p, "inplace dilation l21 is missing")
        self.assertTrue(p[p.index(l12_1)].hasSameOrientation(l12_1),"l12_1 in getBoundary after dilation has incorrect orientation.")
        self.assertTrue(l20 in p, "inplace dilation l20 is missing")
        self.assertTrue(p[p.index(l20)].hasSameOrientation(l20),"l20 in getBoundary after dilation has incorrect orientation.")

        p=s.getBoundaryLoop()
        self.assertTrue(cl1 == p, "inplace dilation s.getBoundaryLoop does not return cl1")
        self.assertTrue(p.hasSameOrientation(cl1),"cl1 in getBoundaryLoop after dilation has incorrect orientation.")

        self.assertTrue(s.getRecombination() == None, "recombination meshing set.")
        self.assertRaises(ValueError,s.setRecombination,-10*DEG)
        s.setRecombination(30*DEG)
        self.assertTrue(s.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertRaises(ValueError,s.setTransfiniteMeshing,orientation="X")

        # now the same but without holes:
        s=PlaneSurface(cl)
        self.assertTrue(s.getTransfiniteMeshing() == None, "transfinite meshing set.")
        self.assertRaises(ValueError,s.setTransfiniteMeshing,orientation="X")
        self.assertRaises(ValueError,s.setTransfiniteMeshing)

        l01.setElementDistribution(3)
        l12_1.setElementDistribution(3)
        l20.setElementDistribution(3)
        s.setTransfiniteMeshing()
        q=s.getTransfiniteMeshing()
        self.assertTrue(not q == None, "transfinite meshing not set.")
        self.assertTrue(q[1] == RuledSurface.LEFT, "orientation is wrong.")
        self.assertTrue(len(q[0])==3, "number of corner points is wrong.")
        self.assertTrue(p0 in q[0], "p0 is missing.")
        self.assertTrue(p1 in q[0], "p1 is missing.")
        self.assertTrue(p2 in q[0], "p2 is missing.")
        s.resetTransfiniteMeshing()
        self.assertTrue(s.getTransfiniteMeshing() == None, "reset transfinite meshing failed.")

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

        self.assertRaises(TypeError,RuledSurface,l01)

        CC0=RuledSurface(cl1)
        s=-CC0

        cl=s.getBoundaryLoop()
        self.assertTrue(cl == cl1, " wrong boundary loops")
        self.assertTrue(cl.hasSameOrientation(-cl1),"cl1 has incorrect orientation.")

        self.assertTrue(s.hasSameOrientation(s),"has not same orientation like itself")
        self.assertTrue(not s.hasSameOrientation(-s),"has same orientation like -itself")

        crvs=s.getBoundary()
        self.assertTrue(len(crvs) == 3, "too many boundary corves.")
        self.assertTrue(l01 in crvs, "l01 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l01)].hasSameOrientation(-l01),"l01 has incorrect orientation.")
        self.assertTrue(l12_1 in crvs, "l21 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l12_1)].hasSameOrientation(-l12_1),"l12_1 has incorrect orientation.")
        self.assertTrue(l20 in crvs, "l20 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l20)].hasSameOrientation(-l20),"l12_1 has incorrect orientation.")

        self.assertTrue(not s.isColocated(p4),"RuledSurface is colocated with point.")
        self.assertTrue(s.isColocated(s),"RuledSurface is not colocated with its self.")
        self.assertTrue(s.isColocated(RuledSurface(cl1)),"RuledSurface is not colocated with its copy.")
        self.assertTrue(not s.isColocated(RuledSurface(cl2)),"RuledSurface is colocated with different length")
        self.assertTrue(not s.isColocated(RuledSurface(cl3)),"RuledSurface is colocated with same length.")

        s.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")

        p=s.getPrimitives()
        self.assertTrue(len(p) == 10, "too many primitives.")
        self.assertTrue(cl1 in p, "cl1 is missing")
        self.assertTrue(l01 in p, "l01 is missing")
        self.assertTrue(l12_1 in p, "l21 is missing")
        self.assertTrue(l20 in p, "l20 is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")

        sp=s.copy()
        self.assertTrue(isinstance(sp,ReverseRuledSurface), "copy returns is not a RuledSurface.")
        self.assertTrue(not sp == s, "copy equals source")
        self.assertTrue(sp.isColocated(s),"copy is not colocated with its source.")
        cbl=sp.getBoundaryLoop()
        self.assertTrue(not cbl == cl1,"copy uses cl1.")
        cp=sp.getPrimitives()
        self.assertTrue(len(cp) == 10, "copy as too many primitives.")
        self.assertTrue(not cl1 in cp, "copy is using cl1")
        self.assertTrue(not l01 in cp, "copy is using l01")
        self.assertTrue(not l12_1 in cp, "copy is using l21")
        self.assertTrue(not l20 in cp, "copy is using l20")
        self.assertTrue(not p0 in cp, "copy is using p0")
        self.assertTrue(not p1 in cp, "copy is using p1")
        self.assertTrue(not p2 in cp, "copy is using p2")
        self.assertTrue(not p3 in cp, "copy is using p3")
        self.assertTrue(not p4 in cp, "copy is using p4")
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
        self.assertTrue(ds.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"dilation is wrong.")
        cbl=ds.getBoundaryLoop()
        self.assertTrue(not cbl == cl1,"dilation uses cl1.")
        cp=ds.getPrimitives()
        self.assertTrue(len(cp) == 10, "dilation as too many primitives.")
        self.assertTrue(not cl1 in cp, "dilation is using cl1")
        self.assertTrue(not l01 in cp, "dilation is using l01")
        self.assertTrue(not l12_1 in cp, "dilation is using l21")
        self.assertTrue(not l20 in cp, "dilation is using l20")
        self.assertTrue(not p0 in cp, "dilation is using p0")
        self.assertTrue(not p1 in cp, "dilation is using p1")
        self.assertTrue(not p2 in cp, "dilation is using p2")
        self.assertTrue(not p3 in cp, "dilation is using p3")
        self.assertTrue(not p4 in cp, "dilation is using p4")

        s.modifyBy(Dilation(-1.))
        self.assertTrue(s.isColocated(RuledSurface(CurveLoop(l01_m,l12_m,l20_m))),"inplace dilation is wrong.")

        p=s.getPrimitives()
        self.assertTrue(len(p) == 10, "inplace dilation has too many primitives.")
        self.assertTrue(cl1 in p, "inplace dilation cl1 is missing")
        self.assertTrue(l01 in p, "inplace dilation l01 is missing")
        self.assertTrue(l12_1 in p, "inplace dilation l21 is missing")
        self.assertTrue(l20 in p, "inplace dilation l20 is missing")
        self.assertTrue(p0 in p, "inplace dilation p0 is missing")
        self.assertTrue(p1 in p, "inplace dilation p1 is missing")
        self.assertTrue(p2 in p, "inplace dilation p2 is missing")
        self.assertTrue(p3 in p, "inplace dilation p3 is missing")
        self.assertTrue(p4 in p, "inplace dilation p4 is missing")

        p=s.getBoundary()
        self.assertTrue(len(p) == 3, "inplace dilation has too many boundary curves.")
        self.assertTrue(l01 in p, "inplace dilation l01 is missing in boundary curves.")
        self.assertTrue(p[p.index(l01)].hasSameOrientation(-l01),"l01 in getBoundary after dilation has incorrect orientation.")
        self.assertTrue(l12_1 in p, "inplace dilation l21 is missing")
        self.assertTrue(p[p.index(l12_1)].hasSameOrientation(-l12_1),"l12_1 in getBoundary after dilation has incorrect orientation.")
        self.assertTrue(l20 in p, "inplace dilation l20 is missing")
        self.assertTrue(p[p.index(l20)].hasSameOrientation(-l20),"l20 in getBoundary after dilation has incorrect orientation.")

        p=s.getBoundaryLoop()
        self.assertTrue(cl1 == p, "inplace dilation s.getBoundaryLoop does not return cl1")
        self.assertTrue(p.hasSameOrientation(-cl1),"cl1 in getBoundaryLoop after dilation has incorrect orientation.")

        self.assertTrue(s.getRecombination() == None, "recombination meshing set.")
        self.assertRaises(ValueError,s.setRecombination,-10*DEG)
        s.setRecombination(30*DEG)
        self.assertTrue(s.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertRaises(ValueError,s.setTransfiniteMeshing,orientation="X")

        # now the same but without holes:
        s=PlaneSurface(cl1)
        self.assertTrue(s.getTransfiniteMeshing() == None, "transfinite meshing set.")
        self.assertRaises(ValueError,s.setTransfiniteMeshing,orientation="X")
        self.assertRaises(ValueError,s.setTransfiniteMeshing)

        l01.setElementDistribution(3)
        l12_1.setElementDistribution(3)
        l20.setElementDistribution(3)
        s.setTransfiniteMeshing()
        q=s.getTransfiniteMeshing()
        self.assertTrue(not q == None, "transfinite meshing not set.")
        self.assertTrue(q[1] == RuledSurface.LEFT, "orientation is wrong.")
        self.assertTrue(len(q[0])==3, "number of corner points is wrong.")
        self.assertTrue(p0 in q[0], "p0 is missing.")
        self.assertTrue(p1 in q[0], "p1 is missing.")
        self.assertTrue(p2 in q[0], "p2 is missing.")
        s.resetTransfiniteMeshing()
        self.assertTrue(s.getTransfiniteMeshing() == None, "reset transfinite meshing failed.")

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

        self.assertRaises(TypeError,PlaneSurface,l4)
        # self.assertRaises(TypeError,PlaneSurface,cl_a,h) activate if check for points not on plane
        # self.assertRaises(TypeError,PlaneSurface,cl,[h_a])

        s=PlaneSurface(cl,holes=[h])

        cl2=s.getBoundaryLoop()
        self.assertTrue(cl == cl2, " wrong boundary loops")
        self.assertTrue(cl.hasSameOrientation(cl2),"cl has incorrect orientation.")

        hs=s.getHoles()
        self.assertTrue(len(hs) == 1, "one holes expected.")
        self.assertTrue(h==hs[0], "h is not defined as hole.")
        self.assertTrue(hs[0].hasSameOrientation(h),"hs has incorrect orientation.")

        self.assertTrue(s.hasSameOrientation(s),"has not same orientation like itself")
        self.assertTrue(not s.hasSameOrientation(-s),"has same orientation like -itself")

        crvs=s.getBoundary()
        self.assertTrue(len(crvs) == 7, "too many boundary corves.")
        self.assertTrue(l0 in crvs, "l0 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l0)].hasSameOrientation(l0),"l0 has incorrect orientation.")
        self.assertTrue(l1 in crvs, "l1 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l1)].hasSameOrientation(l1),"l1 has incorrect orientation.")
        self.assertTrue(l2 in crvs, "l2 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l2)].hasSameOrientation(l2),"l2 has incorrect orientation.")
        self.assertTrue(l3 in crvs, "l3 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l3)].hasSameOrientation(l3),"l3 has incorrect orientation.")
        self.assertTrue(l4 in crvs, "l4 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l4)].hasSameOrientation(l4),"l4 has incorrect orientation.")
        self.assertTrue(l5 in crvs, "l5 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l5)].hasSameOrientation(l5),"l5 has incorrect orientation.")
        self.assertTrue(l6 in crvs, "l6 is missing in boundary")
        self.assertTrue(crvs[crvs.index(l6)].hasSameOrientation(l6),"l6 has incorrect orientation.")

        self.assertTrue(not s.isColocated(p4),"PlaneSurface is colocated with point.")
        self.assertTrue(s.isColocated(s),"PlaneSurface is not colocated with its self.")
        self.assertTrue(s.isColocated(PlaneSurface(cl,holes=[h])),"PlaneSurface is not colocated with its copy.")
        self.assertTrue(not s.isColocated(PlaneSurface(cl)),"PlaneSurface is colocated with PlaneSurface with same boundary but no hole")
        self.assertTrue(not s.isColocated(PlaneSurface(cl_s,holes=[h])),"PlaneSurface is colocated with PlaneSurface with deformed boundary")
        self.assertTrue(not s.isColocated(PlaneSurface(cl,holes=[h2])),"PlaneSurface is colocated with modified hole")

        s.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.assertTrue(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.assertTrue(p6.getLocalScale()==3., "p6 has wrong local scale.")

        p=s.getPrimitives()
        self.assertTrue(len(p) == 17, "too many primitives.")
        self.assertTrue(s in p, "cl is missing")
        self.assertTrue(cl in p, "cl is missing")
        self.assertTrue(h in p, "h is missing")
        self.assertTrue(l0 in p, "l0 is missing")
        self.assertTrue(l1 in p, "l1 is missing")
        self.assertTrue(l2 in p, "l2 is missing")
        self.assertTrue(l3 in p, "l3 is missing")
        self.assertTrue(l4 in p, "l4 is missing")
        self.assertTrue(l5 in p, "l5 is missing")
        self.assertTrue(l6 in p, "l6 is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")
        self.assertTrue(p5 in p, "p5 is missing")
        self.assertTrue(p6 in p, "p6 is missing")

        sp=s.copy()
        self.assertTrue(isinstance(sp,PlaneSurface), "copy returns is not a PlaneSurface.")
        self.assertTrue(not sp == s, "copy equals source")
        self.assertTrue(sp.isColocated(s),"copy is not colocated with its source.")
        cbl=sp.getBoundaryLoop()
        self.assertTrue(not cbl == cl,"copy uses cl1.")
        hs=sp.getHoles()
        self.assertTrue(len(hs)==1,"copy is missing holes.")
        self.assertTrue(not hs[0]== h,"copy contains h as a hole.")
        cp=sp.getPrimitives()
        self.assertTrue(len(cp) == 17, "copy as too many primitives.")
        self.assertTrue(not s in cp, "copy contains s")
        self.assertTrue(not cl in cp, "copy contains cl")
        self.assertTrue(not h in cp, "copy contains h")
        self.assertTrue(not l0 in cp, "copy contains l0")
        self.assertTrue(not l1 in cp, "copy contains l1")
        self.assertTrue(not l2 in cp, "copy contains l2")
        self.assertTrue(not l3 in cp, "copy contains l3")
        self.assertTrue(not l4 in cp, "copy contains l4")
        self.assertTrue(not l5 in cp, "copy contains l5")
        self.assertTrue(not l6 in cp, "copy contains l6")
        self.assertTrue(not p0 in cp, "copy contains p0")
        self.assertTrue(not p1 in cp, "copy contains p1")
        self.assertTrue(not p2 in cp, "copy contains p2")
        self.assertTrue(not p3 in cp, "copy contains p3")
        self.assertTrue(not p4 in cp, "copy contains p4")
        self.assertTrue(not p5 in cp, "copy contains p5")
        self.assertTrue(not p6 in cp, "copy contains p6")
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
        self.assertTrue(ds.isColocated(PlaneSurface(cl_m,holes=[h_m])),"dilation is wrong.")
        cbl=ds.getBoundaryLoop()
        self.assertTrue(not cbl == cl,"dilation uses cl1.")
        hs=ds.getHoles()
        self.assertTrue(len(hs)==1,"dilation is missing holes.")
        self.assertTrue(not hs[0]== h,"dilation contains h as a hole.")
        cp=ds.getPrimitives()
        self.assertTrue(len(cp) == 17, "dilation as too many primitives.")
        self.assertTrue(not s in cp, "dilation contains s")
        self.assertTrue(not cl in cp, "dilation contains cl")
        self.assertTrue(not h in cp, "dilation contains h")
        self.assertTrue(not l0 in cp, "dilation contains l0")
        self.assertTrue(not l1 in cp, "dilation contains l1")
        self.assertTrue(not l2 in cp, "dilation contains l2")
        self.assertTrue(not l3 in cp, "dilation contains l3")
        self.assertTrue(not l4 in cp, "dilation contains l4")
        self.assertTrue(not l5 in cp, "dilation contains l5")
        self.assertTrue(not l6 in cp, "dilation contains l6")
        self.assertTrue(not p0 in cp, "dilation contains p0")
        self.assertTrue(not p1 in cp, "dilation contains p1")
        self.assertTrue(not p2 in cp, "dilation contains p2")
        self.assertTrue(not p3 in cp, "dilation contains p3")
        self.assertTrue(not p4 in cp, "dilation contains p4")
        self.assertTrue(not p5 in cp, "dilation contains p5")
        self.assertTrue(not p6 in cp, "dilation contains p6")

        s.modifyBy(Dilation(-1.))
        self.assertTrue(s.isColocated(PlaneSurface(cl_m,holes=[h_m])),"inplace dilation is wrong.")
        cbl=s.getBoundaryLoop()
        self.assertTrue(cbl == cl,"inplace dilation does not use cl1.")
        self.assertTrue(cl.hasSameOrientation(cbl),"cl has incorrect orientation.")
        hs=s.getHoles()
        self.assertTrue(len(hs)==1,"inplace dilation is missing holes.")
        self.assertTrue(hs[0]== h,"inplace dilation must contain h as a hole.")
        self.assertTrue(hs[0].hasSameOrientation(h),"hole in inplace dilation has incorrect orientation.")

        cp=s.getPrimitives()
        self.assertTrue(len(cp) == 17, "inplace dilation as too many primitives.")
        self.assertTrue(s in cp, "inplace dilation must use s")
        self.assertTrue(cl in cp, "inplace dilation must use cl")
        self.assertTrue(h in cp, "inplace dilation must use h")
        self.assertTrue(l0 in cp, "inplace dilation must use l0")
        self.assertTrue(l1 in cp, "inplace dilation must use l1")
        self.assertTrue(l2 in cp, "inplace dilation must use l2")
        self.assertTrue(l3 in cp, "inplace dilation must use l3")
        self.assertTrue(l4 in cp, "inplace dilation must use l4")
        self.assertTrue(l5 in cp, "inplace dilation must use l5")
        self.assertTrue(l6 in cp, "inplace dilation must use l6")
        self.assertTrue(p0 in cp, "inplace dilation must use p0")
        self.assertTrue(p1 in cp, "inplace dilation must use p1")
        self.assertTrue(p2 in cp, "inplace dilation must use p2")
        self.assertTrue(p3 in cp, "inplace dilation must use p3")
        self.assertTrue(p4 in cp, "inplace dilation must use p4")
        self.assertTrue(p5 in cp, "inplace dilation must use p5")
        self.assertTrue(p6 in cp, "inplace dilation must use p6")

        self.assertTrue(s.getRecombination() == None, "recombination meshing set.")
        self.assertRaises(ValueError,s.setRecombination,-10*DEG)
        s.setRecombination(30*DEG)
        self.assertTrue(s.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertRaises(ValueError,s.setTransfiniteMeshing,orientation="X")

        # now the same but without holes:
        s=PlaneSurface(cl)
        self.assertTrue(s.getTransfiniteMeshing() == None, "transfinite meshing set.")
        self.assertRaises(ValueError,s.setTransfiniteMeshing,orientation="X")
        self.assertRaises(ValueError,s.setTransfiniteMeshing)

        l0.setElementDistribution(6)
        l1.setElementDistribution(3)
        l2.setElementDistribution(6)
        l3.setElementDistribution(3)
        s.setTransfiniteMeshing()
        q=s.getTransfiniteMeshing()
        self.assertTrue(not q == None, "transfinite meshing not set.")
        self.assertTrue(q[1] == RuledSurface.LEFT, "orientation is wrong.")
        self.assertTrue(len(q[0])==4, "number of corner points is wrong.")
        self.assertTrue(p0 in q[0], "p0 is missing.")
        self.assertTrue(p1 in q[0], "p1 is missing.")
        self.assertTrue(p2 in q[0], "p2 is missing.")
        self.assertTrue(p3 in q[0], "p3 is missing.")
        s.resetTransfiniteMeshing()
        self.assertTrue(s.getTransfiniteMeshing() == None, "reset transfinite meshing failed.")

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

        # self.assertRaises(ValueError,SurfaceLoop,s1,s3)
        # self.assertRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5)
        # self.assertRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5,s5)
        s=SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)

        self.assertTrue(s.hasSameOrientation(s),"has not same orientation like itself")
        self.assertTrue(not s.hasSameOrientation(-s),"has same orientation like -itself")

        cc=s.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many curves.")
        self.assertTrue(s1 in cc, "s1 is missing")
        self.assertTrue(s2 in cc, "s2 is missing")
        self.assertTrue(s3 in cc, "s3 is missing")
        self.assertTrue(s4 in cc, "s4 is missing")
        self.assertTrue(s5 in cc, "s5 is missing")
        self.assertTrue(s6 in cc, "s6 is missing")
        self.assertTrue(s7 in cc, "s7 is missing")
        self.assertTrue(s8 in cc, "s8 is missing")
        self.assertTrue(s9 in cc, "s9 is missing")
        self.assertTrue(s10 in cc, "s10 is missing")

        self.assertTrue(not s.isColocated(p4),"SurfaceLoop is colocated with point.")
        self.assertTrue(s.isColocated(s),"SurfaceLoop is not colocated with its self.")
        self.assertTrue(s.isColocated(-s),"SurfaceLoop is not colocated with its reverse.")
        self.assertTrue(s.isColocated(SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)),"SurfaceLoop is not colocated with its copy.")
        self.assertTrue(s.isColocated(SurfaceLoop(-s10, s1,s2,s3,s4,s5,s6,-s7,-s8,-s9)),"SurfaceLoop is not colocated with its copy with shifted points.")
        self.assertTrue(s.isColocated(SurfaceLoop(s1, -s7, s2, -s9, s3, -s8, s4, -s10, s5,s6)),"SurfaceLoop is not colocated with its copy with shuffled points.")
        self.assertTrue(not s.isColocated(SurfaceLoop(s1_v,s2,s3_v,s4,s5,s6)),"SurfaceLoop is colocated with different SurfaceLoop.")

        self.assertTrue(len(s) == 10, "wrong length")

        s.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.assertTrue(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.assertTrue(p6.getLocalScale()==3., "p6 has wrong local scale.")
        self.assertTrue(p7.getLocalScale()==3., "p7 has wrong local scale.")

        p=s.getPrimitives()
        self.assertTrue(len(p) == 63, "too many primitives.")
        self.assertTrue(s in p, "s is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")
        self.assertTrue(p5 in p, "p5 is missing")
        self.assertTrue(p6 in p, "p6 is missing")
        self.assertTrue(p7 in p, "p7 is missing")
        self.assertTrue(q0 in p, "q0 is missing")
        self.assertTrue(q1 in p, "q1 is missing")
        self.assertTrue(q2 in p, "q2 is missing")
        self.assertTrue(q3 in p, "q3 is missing")
        self.assertTrue(q4 in p, "q4 is missing")
        self.assertTrue(q5 in p, "q5 is missing")
        self.assertTrue(q6 in p, "q6 is missing")
        self.assertTrue(q7 in p, "q7 is missing")
        self.assertTrue(l01 in p, "l01 is missing")
        self.assertTrue(l15 in p, "l15 is missing")
        self.assertTrue(l54 in p, "l54 is missing")
        self.assertTrue(l04 in p, "l04 is missing")
        self.assertTrue(l13 in p, "l13 is missing")
        self.assertTrue(l37 in p, "l37 is missing")
        self.assertTrue(l75 in p, "l75 is missing")
        self.assertTrue(l67 in p, "l67 is missing")
        self.assertTrue(l26 in p, "l26 is missing")
        self.assertTrue(l32 in p, "l32 is missing")
        self.assertTrue(l20 in p, "l20 is missing")
        self.assertTrue(l46 in p, "l46 is missing")
        self.assertTrue(m01 in p, "m01 is missing")
        self.assertTrue(m15 in p, "m15 is missing")
        self.assertTrue(m54 in p, "m54 is missing")
        self.assertTrue(m40 in p, "m40 is missing")
        self.assertTrue(m23 in p, "m23 is missing")
        self.assertTrue(m37 in p, "m37 is missing")
        self.assertTrue(m76 in p, "m76 is missing")
        self.assertTrue(m62 in p, "m62 is missing")
        self.assertTrue(m02 in p, "m02 is missing")
        self.assertTrue(m13 in p, "m13 is missing")
        self.assertTrue(m46 in p, "m46 is missing")
        self.assertTrue(m57 in p, "m57 is missing")
        self.assertTrue(cl_l1 in p, "cl_l1 is missing")
        self.assertTrue(cl_m1 in p, "cl_m1 is missing")
        self.assertTrue(s1 in p, "s1 is missing")
        self.assertTrue(cl_l2 in p, "cl_l2 is missing")
        self.assertTrue(s2 in p, "s2 is missing")
        self.assertTrue(cl_l3 in p, "cl_l3 is missing")
        self.assertTrue(cl_m3 in p, "cl_m3 is missing")
        self.assertTrue(s3 in p, "s3 is missing")
        self.assertTrue(cl_l4 in p, "cl_l4 is missing")
        self.assertTrue(s4 in p, "s4 is missing")
        self.assertTrue(cl_l5 in p, "cl_l5 is missing")
        self.assertTrue(s5 in p, "s5 is missing")
        self.assertTrue(cl_l6 in p, "cl_l6 is missing")
        self.assertTrue(s6 in p, "s6 is missing")
        self.assertTrue(cl_m7 in p, "cl_m7 is missing")
        self.assertTrue(s7 in p, "s7 is missing")
        self.assertTrue(cl_m8 in p, "cl_m8 is missing")
        self.assertTrue(s8 in p, "s8 is missing")
        self.assertTrue(cl_m9 in p, "cl_m9 is missing")
        self.assertTrue(s9 in p, "s9 is missing")
        self.assertTrue(cl_m10 in p, "cl_m10 is missing")
        self.assertTrue(s10 in p, "s10 is missing")

        cp=s.copy()
        self.assertTrue(isinstance(cp,SurfaceLoop), "copy returns is not an arc.")
        self.assertTrue(not cp == s, "copy equals source")
        self.assertTrue(cp.isColocated(s),"copy is not colocated with its source.")
        cc=cp.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many primitives in copy.")
        self.assertTrue(not s1 in cc,"copy uses s1.")
        self.assertTrue(not s2 in cc,"copy uses s2.")
        self.assertTrue(not s3 in cc,"copy uses s3.")
        self.assertTrue(not s4 in cc,"copy uses s4.")
        self.assertTrue(not s5 in cc,"copy uses s5.")
        self.assertTrue(not s6 in cc,"copy uses s6.")
        self.assertTrue(not s7 in cc,"copy uses s7.")
        self.assertTrue(not s8 in cc,"copy uses s8.")
        self.assertTrue(not s9 in cc,"copy uses s9.")
        self.assertTrue(not s10 in cc,"copy uses s10.")

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
        self.assertTrue(dc.isColocated(s_m),"dilation is wrong.")
        cc=dc.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many surfaces in copy.")
        self.assertTrue(not s1 in cc,"dilation uses s1.")
        self.assertTrue(not s2 in cc,"dilation uses s2.")
        self.assertTrue(not s3 in cc,"dilation uses s3.")
        self.assertTrue(not s4 in cc,"dilation uses s4.")
        self.assertTrue(not s5 in cc,"dilation uses s5.")
        self.assertTrue(not s6 in cc,"dilation uses s6.")
        self.assertTrue(not s7 in cc,"dilation uses s7.")
        self.assertTrue(not s8 in cc,"dilation uses s8.")
        self.assertTrue(not s9 in cc,"dilation uses s9.")
        self.assertTrue(not s10 in cc,"dilation uses s10.")

        s.modifyBy(Dilation(-1.))
        self.assertTrue(s.isColocated(s_m),"inplace dilation is wrong.")
        cc=s.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many primitives in modified object.")
        self.assertTrue(s1 in cc,"dilation misses s1.")
        self.assertTrue(cc[cc.index(s1)].hasSameOrientation(s1),"s1 in modified object has wrong orientation.")
        self.assertTrue(s1.isColocated(s1_m),"s1 in modified object as wrong location.")
        self.assertTrue(s2 in cc,"dilation misses s2.")
        self.assertTrue(cc[cc.index(s2)].hasSameOrientation(s2),"s2 in modified object has wrong orientation.")
        self.assertTrue(s2.isColocated(s2_m),"s2 in modified object as wrong location.")
        self.assertTrue(s3 in cc,"dilation misses s3.")
        self.assertTrue(cc[cc.index(s3)].hasSameOrientation(s3),"s3 in modified object has wrong orientation.")
        self.assertTrue(s3.isColocated(s3_m),"s3 in modified object as wrong location.")
        self.assertTrue(s4 in cc,"dilation misses s4.")
        self.assertTrue(cc[cc.index(s4)].hasSameOrientation(s4),"s4 in modified object has wrong orientation.")
        self.assertTrue(s4.isColocated(s4_m),"s4 in modified object as wrong location.")
        self.assertTrue(s5 in cc,"dilation misses s5.")
        self.assertTrue(cc[cc.index(s5)].hasSameOrientation(s5),"s5 in modified object has wrong orientation.")
        self.assertTrue(s5.isColocated(s5_m),"s5 in modified object as wrong location.")
        self.assertTrue(s6 in cc,"dilation misses s6.")
        self.assertTrue(cc[cc.index(s6)].hasSameOrientation(s6),"s6 in modified object has wrong orientation.")
        self.assertTrue(s6.isColocated(s6_m),"s6 in modified object as wrong location.")
        self.assertTrue(s7 in cc,"dilation misses s7.")
        self.assertTrue(cc[cc.index(s7)].hasSameOrientation(-s7),"s7 in modified object has wrong orientation.")
        self.assertTrue(s7.isColocated(s7_m),"s7 in modified object as wrong location.")
        self.assertTrue(s8 in cc,"dilation misses s8.")
        self.assertTrue(cc[cc.index(s8)].hasSameOrientation(-s8),"s8 in modified object has wrong orientation.")
        self.assertTrue(s8.isColocated(s8_m),"s8 in modified object as wrong location.")
        self.assertTrue(s9 in cc,"dilation misses s9.")
        self.assertTrue(cc[cc.index(s9)].hasSameOrientation(-s9),"s9 in modified object has wrong orientation.")
        self.assertTrue(s9.isColocated(s9_m),"s9 in modified object as wrong location.")
        self.assertTrue(s10 in cc,"dilation misses s10.")
        self.assertTrue(cc[cc.index(s10)].hasSameOrientation(-s10),"s10 in modified object has wrong orientation.")
        self.assertTrue(s10.isColocated(s10_m),"s10 in modified object as wrong location.")

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

        # self.assertRaises(ValueError,SurfaceLoop,s1,s3)
        # self.assertRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5)
        # self.assertRaises(ValueError,SurfaceLoop,s1,s2,s3,s4,s5,s5)

        CC0=SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)
        s=-CC0

        self.assertTrue(s.hasSameOrientation(s),"has not same orientation like itself")
        self.assertTrue(not s.hasSameOrientation(-s),"has same orientation like -itself")

        cc=s.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many curves.")
        self.assertTrue(s1 in cc, "s1 is missing")
        self.assertTrue(s2 in cc, "s2 is missing")
        self.assertTrue(s3 in cc, "s3 is missing")
        self.assertTrue(s4 in cc, "s4 is missing")
        self.assertTrue(s5 in cc, "s5 is missing")
        self.assertTrue(s6 in cc, "s6 is missing")
        self.assertTrue(s7 in cc, "s7 is missing")
        self.assertTrue(s8 in cc, "s8 is missing")
        self.assertTrue(s9 in cc, "s9 is missing")
        self.assertTrue(s10 in cc, "s10 is missing")

        self.assertTrue(not s.isColocated(p4),"SurfaceLoop is colocated with point.")
        self.assertTrue(s.isColocated(s),"SurfaceLoop is not colocated with its self.")
        self.assertTrue(s.isColocated(-s),"SurfaceLoop is not colocated with its reverse.")
        self.assertTrue(s.isColocated(SurfaceLoop(s1,s2,s3,s4,s5,s6,-s7,-s8,-s9,-s10)),"SurfaceLoop is not colocated with its copy.")
        self.assertTrue(s.isColocated(SurfaceLoop(-s10, s1,s2,s3,s4,s5,s6,-s7,-s8,-s9)),"SurfaceLoop is not colocated with its copy with shifted points.")
        self.assertTrue(s.isColocated(SurfaceLoop(s1, -s7, s2, -s9, s3, -s8, s4, -s10, s5,s6)),"SurfaceLoop is not colocated with its copy with shuffled points.")
        self.assertTrue(not s.isColocated(SurfaceLoop(s1_v,s2,s3_v,s4,s5,s6)),"SurfaceLoop is colocated with different SurfaceLoop.")

        self.assertTrue(len(s) == 10, "wrong length")

        s.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.assertTrue(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.assertTrue(p6.getLocalScale()==3., "p6 has wrong local scale.")
        self.assertTrue(p7.getLocalScale()==3., "p7 has wrong local scale.")

        p=s.getPrimitives()
        self.assertTrue(len(p) == 63, "too many primitives.")
        self.assertTrue(s in p, "s is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")
        self.assertTrue(p5 in p, "p5 is missing")
        self.assertTrue(p6 in p, "p6 is missing")
        self.assertTrue(p7 in p, "p7 is missing")
        self.assertTrue(q0 in p, "q0 is missing")
        self.assertTrue(q1 in p, "q1 is missing")
        self.assertTrue(q2 in p, "q2 is missing")
        self.assertTrue(q3 in p, "q3 is missing")
        self.assertTrue(q4 in p, "q4 is missing")
        self.assertTrue(q5 in p, "q5 is missing")
        self.assertTrue(q6 in p, "q6 is missing")
        self.assertTrue(q7 in p, "q7 is missing")
        self.assertTrue(l01 in p, "l01 is missing")
        self.assertTrue(l15 in p, "l15 is missing")
        self.assertTrue(l54 in p, "l54 is missing")
        self.assertTrue(l04 in p, "l04 is missing")
        self.assertTrue(l13 in p, "l13 is missing")
        self.assertTrue(l37 in p, "l37 is missing")
        self.assertTrue(l75 in p, "l75 is missing")
        self.assertTrue(l67 in p, "l67 is missing")
        self.assertTrue(l26 in p, "l26 is missing")
        self.assertTrue(l32 in p, "l32 is missing")
        self.assertTrue(l20 in p, "l20 is missing")
        self.assertTrue(l46 in p, "l46 is missing")
        self.assertTrue(m01 in p, "m01 is missing")
        self.assertTrue(m15 in p, "m15 is missing")
        self.assertTrue(m54 in p, "m54 is missing")
        self.assertTrue(m40 in p, "m40 is missing")
        self.assertTrue(m23 in p, "m23 is missing")
        self.assertTrue(m37 in p, "m37 is missing")
        self.assertTrue(m76 in p, "m76 is missing")
        self.assertTrue(m62 in p, "m62 is missing")
        self.assertTrue(m02 in p, "m02 is missing")
        self.assertTrue(m13 in p, "m13 is missing")
        self.assertTrue(m46 in p, "m46 is missing")
        self.assertTrue(m57 in p, "m57 is missing")
        self.assertTrue(cl_l1 in p, "cl_l1 is missing")
        self.assertTrue(cl_m1 in p, "cl_m1 is missing")
        self.assertTrue(s1 in p, "s1 is missing")
        self.assertTrue(cl_l2 in p, "cl_l2 is missing")
        self.assertTrue(s2 in p, "s2 is missing")
        self.assertTrue(cl_l3 in p, "cl_l3 is missing")
        self.assertTrue(cl_m3 in p, "cl_m3 is missing")
        self.assertTrue(s3 in p, "s3 is missing")
        self.assertTrue(cl_l4 in p, "cl_l4 is missing")
        self.assertTrue(s4 in p, "s4 is missing")
        self.assertTrue(cl_l5 in p, "cl_l5 is missing")
        self.assertTrue(s5 in p, "s5 is missing")
        self.assertTrue(cl_l6 in p, "cl_l6 is missing")
        self.assertTrue(s6 in p, "s6 is missing")
        self.assertTrue(cl_m7 in p, "cl_m7 is missing")
        self.assertTrue(s7 in p, "s7 is missing")
        self.assertTrue(cl_m8 in p, "cl_m8 is missing")
        self.assertTrue(s8 in p, "s8 is missing")
        self.assertTrue(cl_m9 in p, "cl_m9 is missing")
        self.assertTrue(s9 in p, "s9 is missing")
        self.assertTrue(cl_m10 in p, "cl_m10 is missing")
        self.assertTrue(s10 in p, "s10 is missing")

        cp=s.copy()
        self.assertTrue(isinstance(cp,ReverseSurfaceLoop), "copy returns is not ReverseSurfaceLoop.")
        self.assertTrue(not cp == s, "copy equals source")
        self.assertTrue(cp.isColocated(s),"copy is not colocated with its source.")
        cc=cp.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many primitives in copy.")
        self.assertTrue(not s1 in cc,"copy uses s1.")
        self.assertTrue(not s2 in cc,"copy uses s2.")
        self.assertTrue(not s3 in cc,"copy uses s3.")
        self.assertTrue(not s4 in cc,"copy uses s4.")
        self.assertTrue(not s5 in cc,"copy uses s5.")
        self.assertTrue(not s6 in cc,"copy uses s6.")
        self.assertTrue(not s7 in cc,"copy uses s7.")
        self.assertTrue(not s8 in cc,"copy uses s8.")
        self.assertTrue(not s9 in cc,"copy uses s9.")
        self.assertTrue(not s10 in cc,"copy uses s10.")

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
        self.assertTrue(dc.isColocated(s_m),"dilation is wrong.")
        cc=dc.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many surfaces in copy.")
        self.assertTrue(not s1 in cc,"dilation uses s1.")
        self.assertTrue(not s2 in cc,"dilation uses s2.")
        self.assertTrue(not s3 in cc,"dilation uses s3.")
        self.assertTrue(not s4 in cc,"dilation uses s4.")
        self.assertTrue(not s5 in cc,"dilation uses s5.")
        self.assertTrue(not s6 in cc,"dilation uses s6.")
        self.assertTrue(not s7 in cc,"dilation uses s7.")
        self.assertTrue(not s8 in cc,"dilation uses s8.")
        self.assertTrue(not s9 in cc,"dilation uses s9.")
        self.assertTrue(not s10 in cc,"dilation uses s10.")

        s.modifyBy(Dilation(-1.))
        self.assertTrue(s.isColocated(s_m),"inplace dilation is wrong.")
        cc=s.getSurfaces()
        self.assertTrue(len(cc) == 10, "too many primitives in modified object.")
        self.assertTrue(s1 in cc,"dilation misses s1.")
        self.assertTrue(cc[cc.index(s1)].hasSameOrientation(-s1),"s1 in modified object has wrong orientation.")
        self.assertTrue(s1.isColocated(s1_m),"s1 in modified object as wrong location.")
        self.assertTrue(s2 in cc,"dilation misses s2.")
        self.assertTrue(cc[cc.index(s2)].hasSameOrientation(-s2),"s2 in modified object has wrong orientation.")
        self.assertTrue(s2.isColocated(s2_m),"s2 in modified object as wrong location.")
        self.assertTrue(s3 in cc,"dilation misses s3.")
        self.assertTrue(cc[cc.index(s3)].hasSameOrientation(-s3),"s3 in modified object has wrong orientation.")
        self.assertTrue(s3.isColocated(s3_m),"s3 in modified object as wrong location.")
        self.assertTrue(s4 in cc,"dilation misses s4.")
        self.assertTrue(cc[cc.index(s4)].hasSameOrientation(-s4),"s4 in modified object has wrong orientation.")
        self.assertTrue(s4.isColocated(s4_m),"s4 in modified object as wrong location.")
        self.assertTrue(s5 in cc,"dilation misses s5.")
        self.assertTrue(cc[cc.index(s5)].hasSameOrientation(-s5),"s5 in modified object has wrong orientation.")
        self.assertTrue(s5.isColocated(s5_m),"s5 in modified object as wrong location.")
        self.assertTrue(s6 in cc,"dilation misses s6.")
        self.assertTrue(cc[cc.index(s6)].hasSameOrientation(-s6),"s6 in modified object has wrong orientation.")
        self.assertTrue(s6.isColocated(s6_m),"s6 in modified object as wrong location.")
        self.assertTrue(s7 in cc,"dilation misses s7.")
        self.assertTrue(cc[cc.index(s7)].hasSameOrientation(s7),"s7 in modified object has wrong orientation.")
        self.assertTrue(s7.isColocated(s7_m),"s7 in modified object as wrong location.")
        self.assertTrue(s8 in cc,"dilation misses s8.")
        self.assertTrue(cc[cc.index(s8)].hasSameOrientation(s8),"s8 in modified object has wrong orientation.")
        self.assertTrue(s8.isColocated(s8_m),"s8 in modified object as wrong location.")
        self.assertTrue(s9 in cc,"dilation misses s9.")
        self.assertTrue(cc[cc.index(s9)].hasSameOrientation(s9),"s9 in modified object has wrong orientation.")
        self.assertTrue(s9.isColocated(s9_m),"s9 in modified object as wrong location.")
        self.assertTrue(s10 in cc,"dilation misses s10.")
        self.assertTrue(cc[cc.index(s10)].hasSameOrientation(s10),"s10 in modified object has wrong orientation.")
        self.assertTrue(s10.isColocated(s10_m),"s10 in modified object as wrong location.")

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


        self.assertRaises(TypeError,Volume,s1)
        self.assertRaises(TypeError,Volume,s_out,[s1])
        v=Volume(s_out,holes=[s_inner])

        self.assertRaises(NotImplementedError,v.__neg__)

        cl2=v.getSurfaceLoop()
        self.assertTrue(s_out == cl2, " wrong boundary loops")
        self.assertTrue(s_out.hasSameOrientation(cl2),"cl has incorrect orientation.")

        hs=v.getHoles()
        self.assertTrue(len(hs) == 1, "one holes expected.")
        self.assertTrue(s_inner==hs[0], "h is not defined as hole.")
        self.assertTrue(hs[0].hasSameOrientation(s_inner),"hs has incorrect orientation.")

        cc=v.getBoundary()
        self.assertTrue(len(cc) == 12, "too many curves.")
        self.assertTrue(s1 in cc, "s1 is missing")
        self.assertTrue(s2 in cc, "s2 is missing")
        self.assertTrue(s3 in cc, "s3 is missing")
        self.assertTrue(s4 in cc, "s4 is missing")
        self.assertTrue(s5 in cc, "s5 is missing")
        self.assertTrue(s6 in cc, "s6 is missing")
        self.assertTrue(s1_i in cc, "s1_i is missing")
        self.assertTrue(s2_i in cc, "s2_i is missing")
        self.assertTrue(s3_i in cc, "s3_i is missing")
        self.assertTrue(s4_i in cc, "s4_i is missing")
        self.assertTrue(s5_i in cc, "s5_i is missing")
        self.assertTrue(s6_i in cc, "s6_i is missing")

        self.assertTrue(not v.isColocated(p4),"Volume is colocated with point.")
        self.assertTrue(v.isColocated(v),"Volume is not colocated with its self.")
        self.assertTrue(v.isColocated(Volume(s_out,holes=[s_inner])),"Volume is not colocated with its copy.")
        self.assertTrue(not v.isColocated(Volume(s_out,holes=[0.3*s_inner])),"Volume is not colocated with volume of modifies volume .")

        v.setLocalScale(3.)
        self.assertTrue(p0.getLocalScale()==3., "p0 has wrong local scale.")
        self.assertTrue(p1.getLocalScale()==3., "p1 has wrong local scale.")
        self.assertTrue(p2.getLocalScale()==3., "p2 has wrong local scale.")
        self.assertTrue(p3.getLocalScale()==3., "p3 has wrong local scale.")
        self.assertTrue(p4.getLocalScale()==3., "p4 has wrong local scale.")
        self.assertTrue(p5.getLocalScale()==3., "p5 has wrong local scale.")
        self.assertTrue(p6.getLocalScale()==3., "p6 has wrong local scale.")
        self.assertTrue(p7.getLocalScale()==3., "p7 has wrong local scale.")
        self.assertTrue(p0_i.getLocalScale()==3., "p0_i has wrong local scale.")
        self.assertTrue(p1_i.getLocalScale()==3., "p1_i has wrong local scale.")
        self.assertTrue(p2_i.getLocalScale()==3., "p2_i has wrong local scale.")
        self.assertTrue(p3_i.getLocalScale()==3., "p3_i has wrong local scale.")
        self.assertTrue(p4_i.getLocalScale()==3., "p4_i has wrong local scale.")
        self.assertTrue(p5_i.getLocalScale()==3., "p5_i has wrong local scale.")
        self.assertTrue(p6_i.getLocalScale()==3., "p6_i has wrong local scale.")
        self.assertTrue(p7_i.getLocalScale()==3., "p7_i has wrong local scale.")

        p=v.getPrimitives()
        self.assertTrue(len(p) == 67, "too many primitives.")
        self.assertTrue(v in p, "v is missing")
        self.assertTrue(p0 in p, "p0 is missing")
        self.assertTrue(p1 in p, "p1 is missing")
        self.assertTrue(p2 in p, "p2 is missing")
        self.assertTrue(p3 in p, "p3 is missing")
        self.assertTrue(p4 in p, "p4 is missing")
        self.assertTrue(p5 in p, "p5 is missing")
        self.assertTrue(p6 in p, "p6 is missing")
        self.assertTrue(p7 in p, "p7 is missing")
        self.assertTrue(l01 in p, "l01 is missing")
        self.assertTrue(l15 in p, "l15 is missing")
        self.assertTrue(l54 in p, "l54 is missing")
        self.assertTrue(l40 in p, "l40 is missing")
        self.assertTrue(l23 in p, "l23 is missing")
        self.assertTrue(l37 in p, "l37 is missing")
        self.assertTrue(l76 in p, "l76 is missing")
        self.assertTrue(l62 in p, "l62 is missing")
        self.assertTrue(l13 in p, "l13 is missing")
        self.assertTrue(l57 in p, "l57 is missing")
        self.assertTrue(l02 in p, "l02 is missing")
        self.assertTrue(l46 in p, "l46 is missing")
        self.assertTrue(cl1 in p, "cl1 is missing")
        self.assertTrue(s1 in p, "s1 is missing")
        self.assertTrue(cl2 in p, "cl2 is missing")
        self.assertTrue(s2 in p, "s2 is missing")
        self.assertTrue(cl3 in p, "cl3  is missing")
        self.assertTrue(s3 in p, "s3  is missing")
        self.assertTrue(cl4 in p, "cl4  is missing")
        self.assertTrue(s4 in p, "s4  is missing")
        self.assertTrue(cl5 in p, "cl5  is missing")
        self.assertTrue(s5 in p, "s5  is missing")
        self.assertTrue(cl6 in p, "cl6  is missing")
        self.assertTrue(s6 in p, "s6  is missing")
        self.assertTrue(s_out in p, "s_out is missing")
        self.assertTrue(p0_i in p, "p0_i is missing")
        self.assertTrue(p1_i in p, "p1_i is missing")
        self.assertTrue(p2_i in p, "p2_i is missing")
        self.assertTrue(p3_i in p, "p3_i is missing")
        self.assertTrue(p4_i in p, "p4_i is missing")
        self.assertTrue(p5_i in p, "p5_i is missing")
        self.assertTrue(p6_i in p, "p6_i is missing")
        self.assertTrue(p7_i in p, "p7_i is missing")
        self.assertTrue(l01_i in p, "l01_i is missing")
        self.assertTrue(l15_i in p, "l15_i is missing")
        self.assertTrue(l54_i in p, "l54_i is missing")
        self.assertTrue(l40_i in p, "l40_i is missing")
        self.assertTrue(l23_i in p, "l23_i is missing")
        self.assertTrue(l37_i in p, "l37_i is missing")
        self.assertTrue(l76_i in p, "l76_i is missing")
        self.assertTrue(l62_i in p, "l62_i is missing")
        self.assertTrue(l13_i in p, "l13_i is missing")
        self.assertTrue(l57_i in p, "l57_i is missing")
        self.assertTrue(l02_i in p, "l02_i is missing")
        self.assertTrue(l46_i in p, "l46_i is missing")
        self.assertTrue(cl1_i in p, "cl1_i is missing")
        self.assertTrue(s1_i in p, "s1_i is missing")
        self.assertTrue(cl2_i in p, "cl2_i is missing")
        self.assertTrue(s2_i in p, "s2_i is missing")
        self.assertTrue(cl3_i in p, "cl3_i  is missing")
        self.assertTrue(s3_i in p, "s3_i  is missing")
        self.assertTrue(cl4_i in p, "cl4_i  is missing")
        self.assertTrue(s4_i in p, "s4_i  is missing")
        self.assertTrue(cl5_i in p, "cl5_i  is missing")
        self.assertTrue(s5_i in p, "s5_i  is missing")
        self.assertTrue(cl6_i in p, "cl6_i  is missing")
        self.assertTrue(s6_i in p, "s6_i  is missing")
        self.assertTrue(s_inner in p, "s_inner  is missing")

        cp=v.copy()
        self.assertTrue(isinstance(cp,Volume), "copy returns is not an arc.")
        self.assertTrue(not cp == v, "copy equals source")
        self.assertTrue(cp.isColocated(v),"copy is not colocated with its source.")
        cl2=cp.getSurfaceLoop()
        self.assertTrue(not s_out == cl2, "copy uses s_out")
        hs=cp.getHoles()
        self.assertTrue(len(hs) == 1, "copy: one holes expected.")
        self.assertTrue(not s_inner==hs[0], "copy: uses s_inner.")
        cc=cp.getBoundary()
        self.assertTrue(len(cc) == 12, "too many primitives in copy.")
        self.assertTrue(not s1 in cc,"copy uses s1.")
        self.assertTrue(not s2 in cc,"copy uses s2.")
        self.assertTrue(not s3 in cc,"copy uses s3.")
        self.assertTrue(not s4 in cc,"copy uses s4.")
        self.assertTrue(not s5 in cc,"copy uses s5.")
        self.assertTrue(not s6 in cc,"copy uses s6.")
        self.assertTrue(not s1_i in cc,"copy uses s1_i.")
        self.assertTrue(not s2_i in cc,"copy uses s2_i.")
        self.assertTrue(not s3_i in cc,"copy uses s3_i.")
        self.assertTrue(not s4_i in cc,"copy uses s4_i.")
        self.assertTrue(not s5_i in cc,"copy uses s5_i.")
        self.assertTrue(not s6_i in cc,"copy uses s6_i.")

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
        self.assertTrue(dc.isColocated(v),"dilation is wrong.")
        cl2=cp.getSurfaceLoop()
        self.assertTrue(not s_out == cl2, "copy uses s_out")
        hs=cp.getHoles()
        self.assertTrue(len(hs) == 1, "copy: one holes expected.")
        self.assertTrue(not s_inner==hs[0], "copy: uses s_inner.")
        cc=cp.getBoundary()
        self.assertTrue(len(cc) == 12, "too many primitives in copy.")
        self.assertTrue(not s1 in cc,"copy uses s1.")
        self.assertTrue(not s2 in cc,"copy uses s2.")
        self.assertTrue(not s3 in cc,"copy uses s3.")
        self.assertTrue(not s4 in cc,"copy uses s4.")
        self.assertTrue(not s5 in cc,"copy uses s5.")
        self.assertTrue(not s6 in cc,"copy uses s6.")
        self.assertTrue(not s1_i in cc,"copy uses s1_i.")
        self.assertTrue(not s2_i in cc,"copy uses s2_i.")
        self.assertTrue(not s3_i in cc,"copy uses s3_i.")
        self.assertTrue(not s4_i in cc,"copy uses s4_i.")
        self.assertTrue(not s5_i in cc,"copy uses s5_i.")
        self.assertTrue(not s6_i in cc,"copy uses s6_i.")

        v.modifyBy(Dilation(-1.))
        self.assertTrue(v.isColocated(v_m),"inplace dilation is wrong.")
        cl2=v.getSurfaceLoop()
        self.assertTrue(s_out == cl2, "inplace dilation does not s_out")
        hs=v.getHoles()
        self.assertTrue(len(hs) == 1, "inplace dilation: one holes expected.")
        self.assertTrue(s_inner==hs[0], "inplace dilation does not s_inner as hole")
        cc=v.getBoundary()
        self.assertTrue(len(cc) == 12, "too many primitives in copy.")
        self.assertTrue(s1 in cc,"inplace dilation does not use s1.")
        self.assertTrue(cc[cc.index(s1)].hasSameOrientation(s1),"s1 in modified object has wrong orientation.")
        self.assertTrue(s1.isColocated(s1_m),"s1 in dilated object as wrong location.")
        self.assertTrue(s2 in cc,"inplace dilation does not use s2.")
        self.assertTrue(cc[cc.index(s2)].hasSameOrientation(s2),"s2 in modified object has wrong orientation.")
        self.assertTrue(s2.isColocated(s2_m),"s2 in dilated object as wrong location.")
        self.assertTrue(s3 in cc,"inplace dilation does not use s3.")
        self.assertTrue(cc[cc.index(s3)].hasSameOrientation(s3),"s3 in modified object has wrong orientation.")
        self.assertTrue(s3.isColocated(s3_m),"s3 in dilated object as wrong location.")
        self.assertTrue(s4 in cc,"inplace dilation does not use s4.")
        self.assertTrue(cc[cc.index(s4)].hasSameOrientation(s4),"s4 in modified object has wrong orientation.")
        self.assertTrue(s4.isColocated(s4_m),"s4 in dilated object as wrong location.")
        self.assertTrue(s5 in cc,"inplace dilation does not use s5.")
        self.assertTrue(cc[cc.index(s5)].hasSameOrientation(s5),"s5 in modified object has wrong orientation.")
        self.assertTrue(s5.isColocated(s5_m),"s5 in dilated object as wrong location.")
        self.assertTrue(s6 in cc,"inplace dilation does not use s6.")
        self.assertTrue(cc[cc.index(s6)].hasSameOrientation(s6),"s6 in modified object has wrong orientation.")
        self.assertTrue(s6.isColocated(s6_m),"s6 in dilated object as wrong location.")
        self.assertTrue(s1_i in cc,"inplace dilation does not use s1_i.")
        self.assertTrue(cc[cc.index(s1_i)].hasSameOrientation(s1_i),"s1_i in modified object has wrong orientation.")
        self.assertTrue(s1_i.isColocated(s1_i_m),"s1_i in dilated object as wrong location.")
        self.assertTrue(s2_i in cc,"inplace dilation does not use s2_i.")
        self.assertTrue(cc[cc.index(s2_i)].hasSameOrientation(s2_i),"s2_i in modified object has wrong orientation.")
        self.assertTrue(s2_i.isColocated(s2_i_m),"s2_i in dilated object as wrong location.")
        self.assertTrue(s3_i in cc,"inplace dilation does not use s3_i.")
        self.assertTrue(cc[cc.index(s3_i)].hasSameOrientation(s3_i),"s3_i in modified object has wrong orientation.")
        self.assertTrue(s3_i.isColocated(s3_i_m),"s3_i in dilated object as wrong location.")
        self.assertTrue(s4_i in cc,"inplace dilation does not use s4_i.")
        self.assertTrue(cc[cc.index(s4_i)].hasSameOrientation(s4_i),"s4_i in modified object has wrong orientation.")
        self.assertTrue(s4_i.isColocated(s4_i_m),"s4_i in dilated object as wrong location.")
        self.assertTrue(s5_i in cc,"inplace dilation does not use s5_i.")
        self.assertTrue(cc[cc.index(s5_i)].hasSameOrientation(s5_i),"s5_i in modified object has wrong orientation.")
        self.assertTrue(s5_i.isColocated(s5_i_m),"s5_i in dilated object as wrong location.")
        self.assertTrue(s6_i in cc,"inplace dilation does not use s6_i.")
        self.assertTrue(cc[cc.index(s6_i)].hasSameOrientation(s6_i),"s6_i in modified object has wrong orientation.")
        self.assertTrue(s6_i.isColocated(s6_i_m),"s6_i in dilated object as wrong location.")

        # transfinite meshing
        v=Volume(s_out)
        # l01=Line(p0,p1)
        # l15=Line(p1,p5)
        # l54=Line(p5,p4)
        # l40=Line(p4,p0)
        # l23=Line(p2,p3)
        # l37=Line(p3,p7)
        # l76=Line(p7,p6)
        # l62=Line(p6,p2)
        # l13=Line(p1,p3)
        # l57=Line(p5,p7)
        # l02=Line(p0,p2)
        # l46=Line(p4,p6)
        self.assertRaises(ValueError,v.setRecombination,-10*DEG)
        self.assertTrue(v.getTransfiniteMeshing() == None, "transfinite meshing set.")

        v.setElementDistribution(6)
        self.assertTrue( l01.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l15.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l54.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l40.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l23.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l37.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l76.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l62.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l13.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l57.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l02.getElementDistribution()[0] == 6 , "element distribution wrong")
        self.assertTrue( l46.getElementDistribution()[0] == 6 , "element distribution wrong")

        v.setRecombination(30*DEG)
        self.assertTrue(s1.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertTrue(s2.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertTrue(s3.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertTrue(s4.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertTrue(s5.getRecombination() == 30*DEG, "recombination parameter wrong.")
        self.assertTrue(s6.getRecombination() == 30*DEG, "recombination parameter wrong.")
        # now the same but without holes:
        self.assertRaises(ValueError,v.setTransfiniteMeshing,orientation="X")
        v.setTransfiniteMeshing(RuledSurface.RIGHT)
        q=v.getTransfiniteMeshing()
        self.assertTrue(not q == None, "transfinite meshing not set.")
        self.assertTrue(not s1.getTransfiniteMeshing() == None, "recombination parameter wrong.")
        self.assertTrue(s1.getTransfiniteMeshing()[1] == RuledSurface.RIGHT, "orientation is wrong.")
        self.assertTrue(not s2.getTransfiniteMeshing() == None, "recombination parameter wrong.")
        self.assertTrue(s2.getTransfiniteMeshing()[1] == RuledSurface.RIGHT, "orientation is wrong.")
        self.assertTrue(not s3.getTransfiniteMeshing() == None, "recombination parameter wrong.")
        self.assertTrue(s3.getTransfiniteMeshing()[1] == RuledSurface.RIGHT, "orientation is wrong.")
        self.assertTrue(not s4.getTransfiniteMeshing() == None, "recombination parameter wrong.")
        self.assertTrue(s4.getTransfiniteMeshing()[1] == RuledSurface.RIGHT, "orientation is wrong.")
        self.assertTrue(not s5.getTransfiniteMeshing() == None, "recombination parameter wrong.")
        self.assertTrue(s5.getTransfiniteMeshing()[1] == RuledSurface.RIGHT, "orientation is wrong.")
        self.assertTrue(not s6.getTransfiniteMeshing() == None, "recombination parameter wrong.")
        self.assertTrue(s6.getTransfiniteMeshing()[1] == RuledSurface.RIGHT, "orientation is wrong.")

   def test_PropertySet1D(self):
       p0=Point(1.,2.,3.,local_scale=9.)
       p1=Point(0.,0.,0.,local_scale=9.)
       p3=Point(8.,6,6,local_scale=9.)
       p4=Point(8.,6,-6,local_scale=9.)

       l0=Line(p0,p1)
       l1=Arc(p3,p1,p4)
       l2=Line(p4,p0)
       # create property set with dim:
       self.assertRaises(TypeError,PropertySet,"test0",l0, p0)

       ps=PropertySet("test0", l0, l1)
       self.assertTrue(ps.getManifoldClass() == Manifold1D, "wrong manifold")
       self.assertTrue(ps.getDim() == 1, "wrong dimension")
       self.assertTrue(ps.getName() == "test0", "wrong name")
       ps.setName("test1")
       self.assertTrue(ps.getName() == "test1", "wrong new name")

       self.assertTrue(ps.getTag() == 9, "wrong tag")

       self.assertRaises(TypeError,ps.addItem, p0)
       ps.addItem(l2)
       pp=ps.getItems()
       self.assertTrue(len(pp) == 3, "wrong number of items")
       self.assertTrue(l0 in pp, "l0 missing in items.")
       self.assertTrue(l1 in pp, "l1 missing in items.")
       self.assertTrue(l2 in pp, "l2 missing in items.")

       pp=ps.getPrimitives()
       self.assertTrue(len(pp) == 8, "wrong number of items")
       self.assertTrue(ps in pp, "ps missing in items.")
       self.assertTrue(l0 in pp, "l0 missing in items.")
       self.assertTrue(l1 in pp, "l1 missing in items.")
       self.assertTrue(l2 in pp, "l2 missing in items.")
       self.assertTrue(p0 in pp, "p0 missing in items.")
       self.assertTrue(p1 in pp, "p1 missing in items.")
       self.assertTrue(p3 in pp, "p3 missing in items.")
       self.assertTrue(p4 in pp, "p4 missing in items.")

       ps.clearItems()
       self.assertTrue(len(ps.getItems()) == 0, "cleaning items failed.")

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
       self.assertRaises(TypeError,PropertySet,"test0",s, p0)

       ps=PropertySet("test0", s)
       self.assertTrue(ps.getManifoldClass() == Manifold2D, "wrong manifold")
       self.assertTrue(ps.getDim() == 2, "wrong dimension")
       self.assertTrue(ps.getName() == "test0", "wrong name")
       ps.setName("test1")
       self.assertTrue(ps.getName() == "test1", "wrong new name")

       self.assertTrue(ps.getTag() == 19, "wrong tag")

       pp=ps.getPrimitives()
       self.assertTrue(len(pp) == 18, "wrong number of items")
       self.assertTrue(ps in pp, "ps missing in items.")
       self.assertTrue(s in pp, "s missing in items.")
       self.assertTrue(h in pp, "h missing in items.")
       self.assertTrue(cl in pp, "cl missing in items.")
       self.assertTrue(l0 in pp, "l0 missing in items.")
       self.assertTrue(l1 in pp, "l1 missing in items.")
       self.assertTrue(l2 in pp, "l2 missing in items.")
       self.assertTrue(l3 in pp, "l3 missing in items.")
       self.assertTrue(l4 in pp, "l4 missing in items.")
       self.assertTrue(l5 in pp, "l5 missing in items.")
       self.assertTrue(l6 in pp, "l6 missing in items.")
       self.assertTrue(p0 in pp, "p0 missing in items.")
       self.assertTrue(p1 in pp, "p1 missing in items.")
       self.assertTrue(p3 in pp, "p3 missing in items.")
       self.assertTrue(p4 in pp, "p4 missing in items.")
       self.assertTrue(p5 in pp, "p5 missing in items.")
       self.assertTrue(p6 in pp, "p6 missing in items.")

       self.assertRaises(TypeError,ps.addItem, p0)
       pp=ps.getItems()
       self.assertTrue(len(pp) == 1, "wrong number of items")
       self.assertTrue(s in pp, "s missing in items.")

       ps.clearItems()
       self.assertTrue(len(ps.getItems()) == 0, "cleaning items failed.")
       ps.addItem(s)
       self.assertTrue(len(ps.getItems()) == 1, "missing added item")
       self.assertTrue(s in ps.getItems(), "missing added item")

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
       self.assertRaises(TypeError,PropertySet,"test0",v, p0)

       ps=PropertySet("test0", v)
       self.assertTrue(ps.getManifoldClass() == Manifold3D, "wrong manifold")
       self.assertTrue(ps.getDim() == 3, "wrong dimension")
       self.assertTrue(ps.getName() == "test0", "wrong name")
       ps.setName("test1")
       self.assertTrue(ps.getName() == "test1", "wrong new name")

       self.assertTrue(ps.getTag() == 69, "wrong tag")

       pp=ps.getPrimitives()
       self.assertTrue(len(pp) == 68, "too many primitives.")
       self.assertTrue(ps in pp, "ps is missing")
       self.assertTrue(v in pp, "v is missing")
       self.assertTrue(p0 in pp, "p0 is missing")
       self.assertTrue(p1 in pp, "p1 is missing")
       self.assertTrue(p2 in pp, "p2 is missing")
       self.assertTrue(p3 in pp, "p3 is missing")
       self.assertTrue(p4 in pp, "p4 is missing")
       self.assertTrue(p5 in pp, "p5 is missing")
       self.assertTrue(p6 in pp, "p6 is missing")
       self.assertTrue(p7 in pp, "p7 is missing")
       self.assertTrue(l01 in pp, "l01 is missing")
       self.assertTrue(l15 in pp, "l15 is missing")
       self.assertTrue(l54 in pp, "l54 is missing")
       self.assertTrue(l40 in pp, "l40 is missing")
       self.assertTrue(l23 in pp, "l23 is missing")
       self.assertTrue(l37 in pp, "l37 is missing")
       self.assertTrue(l76 in pp, "l76 is missing")
       self.assertTrue(l62 in pp, "l62 is missing")
       self.assertTrue(l13 in pp, "l13 is missing")
       self.assertTrue(l57 in pp, "l57 is missing")
       self.assertTrue(l02 in pp, "l02 is missing")
       self.assertTrue(l46 in pp, "l46 is missing")
       self.assertTrue(cl1 in pp, "cl1 is missing")
       self.assertTrue(s1 in pp, "s1 is missing")
       self.assertTrue(cl2 in pp, "cl2 is missing")
       self.assertTrue(s2 in pp, "s2 is missing")
       self.assertTrue(cl3 in pp, "cl3  is missing")
       self.assertTrue(s3 in pp, "s3  is missing")
       self.assertTrue(cl4 in pp, "cl4  is missing")
       self.assertTrue(s4 in pp, "s4  is missing")
       self.assertTrue(cl5 in pp, "cl5  is missing")
       self.assertTrue(s5 in pp, "s5  is missing")
       self.assertTrue(cl6 in pp, "cl6  is missing")
       self.assertTrue(s6 in pp, "s6  is missing")
       self.assertTrue(s_out in pp, "s_out is missing")
       self.assertTrue(p0_i in pp, "p0_i is missing")
       self.assertTrue(p1_i in pp, "p1_i is missing")
       self.assertTrue(p2_i in pp, "p2_i is missing")
       self.assertTrue(p3_i in pp, "p3_i is missing")
       self.assertTrue(p4_i in pp, "p4_i is missing")
       self.assertTrue(p5_i in pp, "p5_i is missing")
       self.assertTrue(p6_i in pp, "p6_i is missing")
       self.assertTrue(p7_i in pp, "p7_i is missing")
       self.assertTrue(l01_i in pp, "l01_i is missing")
       self.assertTrue(l15_i in pp, "l15_i is missing")
       self.assertTrue(l54_i in pp, "l54_i is missing")
       self.assertTrue(l40_i in pp, "l40_i is missing")
       self.assertTrue(l23_i in pp, "l23_i is missing")
       self.assertTrue(l37_i in pp, "l37_i is missing")
       self.assertTrue(l76_i in pp, "l76_i is missing")
       self.assertTrue(l62_i in pp, "l62_i is missing")
       self.assertTrue(l13_i in pp, "l13_i is missing")
       self.assertTrue(l57_i in pp, "l57_i is missing")
       self.assertTrue(l02_i in pp, "l02_i is missing")
       self.assertTrue(l46_i in pp, "l46_i is missing")
       self.assertTrue(cl1_i in pp, "cl1_i is missing")
       self.assertTrue(s1_i in pp, "s1_i is missing")
       self.assertTrue(cl2_i in pp, "cl2_i is missing")
       self.assertTrue(s2_i in pp, "s2_i is missing")
       self.assertTrue(cl3_i in pp, "cl3_i  is missing")
       self.assertTrue(s3_i in pp, "s3_i  is missing")
       self.assertTrue(cl4_i in pp, "cl4_i  is missing")
       self.assertTrue(s4_i in pp, "s4_i  is missing")
       self.assertTrue(cl5_i in pp, "cl5_i  is missing")
       self.assertTrue(s5_i in pp, "s5_i  is missing")
       self.assertTrue(cl6_i in pp, "cl6_i  is missing")
       self.assertTrue(s6_i in pp, "s6_i  is missing")
       self.assertTrue(s_inner in pp, "s_inner  is missing")

       self.assertRaises(TypeError,ps.addItem, p0)
       pp=ps.getItems()
       self.assertTrue(len(pp) == 1, "wrong number of items")
       self.assertTrue(v in pp, "s missing in items.")

       ps.clearItems()
       self.assertTrue(len(ps.getItems()) == 0, "cleaning items failed.")
       ps.addItem(v)
       self.assertTrue(len(ps.getItems()) == 1, "missing added item")
       self.assertTrue(v in ps.getItems(), "missing added item")

class Test_PyCAD_Design(unittest.TestCase):
   def setUp(self):
         resetGlobalPrimitiveIdCounter()

   def test_tagMap(self):
       self.assertRaises(TypeError,TagMap, { "x" : 5 } )
       self.assertRaises(TypeError,TagMap, { 5 : 10 } )

       m=TagMap({ 5 : "x" })

       m.setMap(x=4,a=6)
       m.setMap(b=6,c=1)
       t=m.getTags()
       self.assertTrue(len(t) == 4)
       self.assertTrue(1 in t)
       self.assertTrue(6 in t)
       self.assertTrue(5 in t)
       self.assertTrue(4 in t)
       self.assertTrue(m.getTags("c") == [1])
       self.assertTrue(m.getTags("x") == [4, 5])
       self.assertTrue(m.getTags("b") == [6])

       self.assertTrue(m.getName(1) == "c")
       self.assertTrue(m.getName(4) == "x")
       self.assertTrue(m.getName(5) == "x")
       self.assertTrue(m.getName(6) == "b")

       t=m.getMapping()
       self.assertTrue(len(t)==4)
       self.assertTrue(1 in t)
       self.assertTrue(4 in t)
       self.assertTrue(5 in t)
       self.assertTrue(6 in t)
       self.assertTrue(t[1] == "c")
       self.assertTrue(t[4] == "x")
       self.assertTrue(t[5] == "x")
       self.assertTrue(t[6] == "b")

       d=m.map(c=10, x = -10, b =60, default =1)
       self.assertTrue(d[1] == 10)
       self.assertTrue(d[4] == -10)
       self.assertTrue(d[5] == -10)
       self.assertTrue(d[6] == 60)

       d=m.map(c=10, x = -10, default =60)
       self.assertTrue(d[1] == 10)
       self.assertTrue(d[4] == -10)
       self.assertTrue(d[5] == -10)
       self.assertTrue(d[6] == 60)

       mxml=m.writeXML()
       m2=TagMap()
       m2.fillFromXML(mxml)
       t=m2.getMapping()
       self.assertTrue(len(t)==4)
       self.assertTrue(1 in t)
       self.assertTrue(4 in t)
       self.assertTrue(5 in t)
       self.assertTrue(6 in t)
       self.assertTrue(t[1] == "c")
       self.assertTrue(t[4] == "x")
       self.assertTrue(t[5] == "x")
       self.assertTrue(t[6] == "b")

   def test_Design(self):
       d=AbstractDesign(dim=2, element_size=0.01, order=1, keep_files=False)
       # check dimension:
       self.assertRaises(ValueError,d.setDim,4)
       d.setDim(3)
       self.assertTrue(d.getDim() == 3)
       # check element order
       self.assertRaises(ValueError,d.setElementOrder,4)
       d.setElementOrder(2)
       self.assertTrue(d.getElementOrder() == 2)
       # check element size
       self.assertRaises(ValueError,d.setElementSize,0)
       d.setElementSize(0.02)
       self.assertTrue(d.getElementSize() == 0.02)
       # files:
       d.setKeepFilesOff()
       self.assertTrue(not d.keepFiles())
       d.setKeepFilesOn()
       self.assertTrue(d.keepFiles())
       # mesh handler:
       self.assertRaises(NotImplementedError,d.getMeshHandler)

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
       self.assertRaises(TypeError,d.addItems,1.)
       d.addItems(ps, pl1, pl2)
       i=d.getItems()
       self.assertTrue(isinstance(i,list))
       self.assertTrue(len(i)==3)
       self.assertTrue(ps in i)
       self.assertTrue(pl1 in i)
       self.assertTrue(pl2 in i)

       p=d.getAllPrimitives()
       self.assertTrue(isinstance(p,list))
       self.assertTrue(len(p)==13)
       self.assertTrue(p0 in p)
       self.assertTrue(p1 in p)
       self.assertTrue(p2 in p)
       self.assertTrue(p3 in p)
       self.assertTrue(l01 in p)
       self.assertTrue(l12 in p)
       self.assertTrue(l23 in p)
       self.assertTrue(l30 in p)
       self.assertTrue(c in p)
       self.assertTrue(s in p)
       self.assertTrue(ps in p)
       self.assertTrue(pl1 in p)
       self.assertTrue(pl2 in p)
       # get tag maps:
       m1=d.getTagMap()
       mm1=m1.getMapping()
       self.assertTrue(len(mm1)==3)
       self.assertTrue(mm1[12]=="A")
       self.assertTrue(mm1[13]=="B")
       self.assertTrue(mm1[11]=="XXXX")
       # clear things:
       d.clearItems()
       i=d.getItems()
       self.assertTrue(isinstance(i,list))
       self.assertTrue(len(i)==0)

   def test_GMSH(self):
       d=GMSHDesign(dim=2, element_size=0.01, order=1, keep_files=False)

       script_name=d.getScriptFileName()
       self.assertTrue(isinstance(script_name,str))
       self.assertTrue(script_name.split(".")[-1] == "geo")
       script_name=os.path.join(PYCAD_WORKDIR,"script.geo")
       d.setScriptFileName(script_name)
       self.assertTrue(script_name == d.getScriptFileName())

       mesh_name=d.getMeshFileName()
       self.assertTrue(isinstance(mesh_name,str))
       self.assertTrue(mesh_name.split(".")[-1] == "msh")
       mesh_name=os.path.join(PYCAD_WORKDIR,"mesh.msh")
       d.setMeshFileName(mesh_name)
       self.assertTrue(mesh_name == d.getMeshFileName())

       d.setOptions(algorithm=d.TETGEN,optimize_quality=False,smoothing=4)
       #cmd=d.getCommandString()
       #self.assertTrue("gmsh -format msh -2 -order 1 -v 3 -o '%s' '%s'"%(os.path.join(".","mesh.msh"), os.path.join(".","script.geo")) == cmd%os.path.join(".","script.geo"))

       d.setOptions(optimize_quality=True)
       #cmd=d.getCommandString()
       #self.assertTrue("gmsh -format msh -2 -order 1 -v 3 -o '%s' '%s'"%(os.path.join(".","mesh.msh"), os.path.join(".","script.geo")) == cmd%os.path.join(".","script.geo"))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e-02};
Point(2) = {1.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e-02};
Point(3) = {1.00000000000000e+00, 1.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e-02};
Point(4) = {0.00000000000000e+00, 1.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e-02};
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};
Physical Surface(11) = {10};
Physical Line(12) = {5, 8};
Physical Line(13) = {6, 7};
Physical Surface(14) = {10};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def not_yet_test_Triangle(self):
       d=TriangleDesign(dim=2, keep_files=False)

       script_name=d.getScriptFileName()
       self.assertTrue(isinstance(script_name,str))
       self.assertTrue(script_name.split(".")[-1] == "poly")
       script_name=os.path.join(PYCAD_WORKDIR,"script.poly")
       d.setScriptFileName(script_name)
       self.assertTrue(script_name == d.getScriptFileName())

       mesh_name=d.getMeshFileName()
       self.assertTrue(isinstance(mesh_name,str))
       mesh_name=os.path.join(PYCAD_WORKDIR,"mesh")
       d.setMeshFileName(mesh_name)
       self.assertTrue(mesh_name == d.getMeshFileName())

       d.setOptions(cmdLineArgs="-Qpqa7.5")
       cmd=d.getCommandString()
       self.assertTrue("triangle -Qpqa7.5 .%s"%(os.path.join(".","script.poly")) == cmd)

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
       self.assertEqual(scrpt, ref)

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(3) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 3.00000000000000e-03};
Point(4) = {3.00000000000000e+00, 3.00000000000000e+00, 3.00000000000000e+00 , 4.00000000000000e-03};
Spline(5) = {1, 2, 3, 4};
Physical Line(6) = {5};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(3) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 3.00000000000000e-03};
Point(4) = {3.00000000000000e+00, 3.00000000000000e+00, 3.00000000000000e+00 , 4.00000000000000e-03};
Spline(5) = {1, 2, 3, 4};
Physical Line(6) = {5};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(3) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 3.00000000000000e-03};
Point(4) = {3.00000000000000e+00, 3.00000000000000e+00, 3.00000000000000e+00 , 4.00000000000000e-03};
Bezier(5) = {1, 2, 3, 4};
Physical Line(6) = {5};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_BSpline(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       p0=Point(0,0,0,0.1)
       p1=Point(1,1,1,0.2)
       p2=Point(2,2,2,0.3)
       p3=Point(3,3,3,0.4)
       self.assertRaises(ValueError,BSpline,p0)
       d.addItems(BSpline(p0,p1,p2,p3))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00, 2.00000000000000e-03};
Point(3) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e-03};
Point(4) = {3.00000000000000e+00, 3.00000000000000e+00, 3.00000000000000e+00, 4.00000000000000e-03};
BSpline(5) = {1, 2, 3, 4};
Physical Line(6) = {5};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00, 2.00000000000000e-03};
Point(3) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e-03};
Point(4) = {3.00000000000000e+00, 3.00000000000000e+00, 3.00000000000000e+00, 4.00000000000000e-03};
BSpline(6) = {1, 2, 3, 4};
Physical Line(7) = {6};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_LineSegment(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       p0=Point(0,0,0,0.1)
       p1=Point(1,1,1,0.2)
       d.addItems(Line(p0,p1))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03 };
Point(2) = {1.00000000000000e+00 , 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03 };
Line(3) = {1, 2};
Physical Line(4) = {3};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_ReverseLineSegment(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       p0=Point(0,0,0,0.1)
       p1=Point(1,1,1,0.2)
       CC0=Line(p0,p1)
       d.addItems(-CC0)

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Line(3) = {1, 2};
Physical Line(4) = {3};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_Arc(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       center=Point(0,0,0,0.1)
       p_start=Point(1,1,1,0.2)
       p_end=Point(1,2,3)
       d.addItems(Arc(center,p_start,p_end))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00, 2.00000000000000e-03};
Point(3) = {1.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e+00, 1.00000000000000e-02};
Circle(4) = {2, 1, 3};
Physical Line(5) = {4};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(3) = {1.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e+00 , 1.00000000000000e-02};
Circle(4) = {2, 1, 3};
Physical Line(5) = {4};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_Ellipse(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       center=Point(0,0,0,0.1)
       mainax=Point(0,1,0,0.1)
       p_start=Point(1,1,1,0.2)
       p_end=Point(1,2,3)
       d.addItems(Ellipse(center,mainax,p_start,p_end))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {0.00000000000000e+00, 1.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(3) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(4) = {1.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e+00 , 1.00000000000000e-02};
Ellipse(5) = {3, 1, 2, 4};
Physical Line(6) = {5};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_ReverseEllipse(self):
       d=GMSHDesign(dim=2, element_size=0.01)
       center=Point(0,0,0,0.1)
       mainax=Point(0,1,0,0.1)
       p_start=Point(1,1,1,0.2)
       p_end=Point(1,2,3)
       CC0=Ellipse(center,mainax,p_start,p_end)
       d.addItems(-CC0)

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {0.00000000000000e+00, 1.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(3) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(4) = {1.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e+00 , 1.00000000000000e-02};
Ellipse(5) = {3, 1, 2, 4};
Physical Line(6) = {5};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(3) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 3.00000000000000e-03};
Point(4) = {3.00000000000000e+00, 3.00000000000000e+00, 3.00000000000000e+00 , 4.00000000000000e-03};
Point(5) = {1.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e+00 , 1.00000000000000e-02};
Line(7) = {1, 2};
Circle(8) = {2, 4, 3};
Spline(12) = {3, 5, 1};
Line Loop(13) = {7, 8, 12};
Ruled Surface(16) = {13};
Physical Surface(17) = {16};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {1.00000000000000e+00, 1.00000000000000e+00, 1.00000000000000e+00 , 2.00000000000000e-03};
Point(3) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 3.00000000000000e-03};
Point(4) = {3.00000000000000e+00, 3.00000000000000e+00, 3.00000000000000e+00 , 4.00000000000000e-03};
Point(5) = {1.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e+00 , 1.00000000000000e-02};
Line(7) = {1, 2};
Circle(8) = {2, 4, 3};
Spline(12) = {3, 5, 1};
Line Loop(13) = {7, 8, 12};
Ruled Surface(16) = {13};
Physical Surface(17) = {16};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_PlaneSurface(self):
       d=GMSHDesign(dim=2, element_size=0.1)
       p0=Point(0,0,0,0.1)
       p1=Point(10,0,0,0.2)
       p2=Point(10,10,0,0.3)
       p3=Point(0,10,3,0.4)
       p4=Point(5,5,0,0.01)
       p5=Point(7,5,0,0.01)
       p6=Point(5,7,0,0.01)
       p7=Point(8,8,0,0.01)
       p8=Point(9,9,0,0.01)

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-02};
Point(2) = {1.00000000000000e+01, 0.00000000000000e+00, 0.00000000000000e+00 , 2.00000000000000e-02};
Point(3) = {1.00000000000000e+01, 1.00000000000000e+01, 0.00000000000000e+00 , 3.00000000000000e-02};
Point(4) = {0.00000000000000e+00, 1.00000000000000e+01, 3.00000000000000e+00 , 4.00000000000000e-02};
Point(5) = {5.00000000000000e+00, 5.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(6) = {7.00000000000000e+00, 5.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
Point(7) = {5.00000000000000e+00, 7.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-03};
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
Physical Surface(30) = {29};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_PlaneSurfaceTransfinite(self):
       d=GMSHDesign(dim=2, element_size=0.1)
       p0=Point(0,0,0,0.1)
       p1=Point(10,0,0,0.2)
       p2=Point(10,10,0,0.3)
       p3=Point(0,10,0,0.4)
       l0=Line(p0,p1)
       l1=Line(p1,p2)
       l2=Line(p2,p3)
       l3=Line(p3,p0)
       l0.setElementDistribution(4)
       l1.setElementDistribution(5)
       l2.setElementDistribution(4)
       l3.setElementDistribution(5)
       p=PlaneSurface(CurveLoop(l0,l1,l2,l3))
       p.setRecombination(30*DEG)
       p.setTransfiniteMeshing()
       d.addItems(p)
       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-02};
Point(2) = {1.00000000000000e+01, 0.00000000000000e+00, 0.00000000000000e+00 , 2.00000000000000e-02};
Point(3) = {1.00000000000000e+01, 1.00000000000000e+01, 0.00000000000000e+00 , 3.00000000000000e-02};
Point(4) = {0.00000000000000e+00, 1.00000000000000e+01, 0.00000000000000e+00 , 4.00000000000000e-02};
Line(5) = {1, 2};
Transfinite Line{5} = 4 Using Progression 1;
Line(6) = {2, 3};
Transfinite Line{6} = 5 Using Progression 1;
Line(7) = {3, 4};
Transfinite Line{7} = 4 Using Progression 1;
Line(8) = {4, 1};
Transfinite Line{8} = 5 Using Progression 1;
Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};
Transfinite Surface{10} = {4,1,2,3} Left;
Recombine Surface {10} = 30.000000;
Physical Surface(11) = {10};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))


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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {-2.00000000000000e+00, -2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {2.00000000000000e+00, -2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(3) = {-2.00000000000000e+00, 2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(4) = {2.00000000000000e+00, 2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(5) = {-2.00000000000000e+00, -2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
Point(6) = {2.00000000000000e+00, -2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
Point(7) = {-2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
Point(8) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
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
Point(34) = {-1.00000000000000e+00 , -1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(35) = {1.00000000000000e+00 , -1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(36) = {-1.00000000000000e+00 , 1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(37) = {1.00000000000000e+00 , 1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(38) = {-1.00000000000000e+00 , -1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
Point(39) = {1.00000000000000e+00 , -1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
Point(40) = {-1.00000000000000e+00 , 1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
Point(41) = {1.00000000000000e+00 , 1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
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
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_VolumeTransfinite(self):
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

       v=Volume(s_out)
       v.setElementDistribution(5)
       v.setRecombination()
       v.setTransfiniteMeshing()
       d.addItems(v)

       scrpt=d.getScriptString();
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {-2.00000000000000e+00, -2.00000000000000e+00, -2.00000000000000e+00, 1.00000000000000e-03};
Point(2) = {2.00000000000000e+00, -2.00000000000000e+00, -2.00000000000000e+00, 1.00000000000000e-03};
Point(3) = {-2.00000000000000e+00, 2.00000000000000e+00, -2.00000000000000e+00, 1.00000000000000e-03};
Point(4) = {2.00000000000000e+00, 2.00000000000000e+00, -2.00000000000000e+00, 1.00000000000000e-03};
Point(5) = {-2.00000000000000e+00, -2.00000000000000e+00, 2.00000000000000e+00, 1.00000000000000e-03};
Point(6) = {2.00000000000000e+00, -2.00000000000000e+00, 2.00000000000000e+00, 1.00000000000000e-03};
Point(7) = {-2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00, 1.00000000000000e-03};
Point(8) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00, 1.00000000000000e-03};
Line(9) = {1, 2};
Transfinite Line{9} = 5 Using Progression 1;
Line(10) = {2, 6};
Transfinite Line{10} = 5 Using Progression 1;
Line(11) = {6, 5};
Transfinite Line{11} = 5 Using Progression 1;
Line(12) = {5, 1};
Transfinite Line{12} = 5 Using Progression 1;
Line(13) = {3, 4};
Transfinite Line{13} = 5 Using Progression 1;
Line(14) = {4, 8};
Transfinite Line{14} = 5 Using Progression 1;
Line(15) = {8, 7};
Transfinite Line{15} = 5 Using Progression 1;
Line(16) = {7, 3};
Transfinite Line{16} = 5 Using Progression 1;
Line(17) = {2, 4};
Transfinite Line{17} = 5 Using Progression 1;
Line(18) = {6, 8};
Transfinite Line{18} = 5 Using Progression 1;
Line(19) = {1, 3};
Transfinite Line{19} = 5 Using Progression 1;
Line(20) = {5, 7};
Transfinite Line{20} = 5 Using Progression 1;
Line Loop(21) = {9, 10, 11, 12};
Plane Surface(22) = {21};
Transfinite Surface{22} = {5,1,2,6} Left;
Recombine Surface {22} = 45.000000;
Line Loop(23) = {13, 14, 15, 16};
Plane Surface(24) = {-23};
Transfinite Surface{24} = {8,4,3,7} Left;
Recombine Surface {24} = 45.000000;
Line Loop(25) = {17, 14, -18, -10};
Plane Surface(26) = {25};
Transfinite Surface{26} = {6,2,4,8} Left;
Recombine Surface {26} = 45.000000;
Line Loop(27) = {20, 16, -19, -12};
Plane Surface(28) = {-27};
Transfinite Surface{28} = {3,7,5,1} Left;
Recombine Surface {28} = 45.000000;
Line Loop(29) = {-9, 19, 13, -17};
Plane Surface(30) = {-29};
Transfinite Surface{30} = {3,1,2,4} Left;
Recombine Surface {30} = 45.000000;
Line Loop(31) = {-11, 18, 15, -20};
Plane Surface(32) = {-31};
Transfinite Surface{32} = {8,6,5,7} Left;
Recombine Surface {32} = 45.000000;
Surface Loop(33) = {22, 24, 26, 28, 30, 32};
Volume(34) = {33};
Transfinite Volume{34};
Physical Volume(35) = {34};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {1.00000000000000e+00, 2.00000000000000e+00, 3.00000000000000e+00 , 9.00000000000000e-02};
Point(2) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 9.00000000000000e-02};
Point(3) = {8.00000000000000e+00, 6.00000000000000e+00, 6.00000000000000e+00 , 9.00000000000000e-02};
Point(4) = {8.00000000000000e+00, 6.00000000000000e+00, -6.00000000000000e+00 , 9.00000000000000e-02};
Line(5) = {1, 2};
Circle(6) = {2, 3, 4};
Line(7) = {4, 1};
Physical Line(8) = {5, 6, 7};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_generate_PropertySet2D(self):
       d=GMSHDesign(dim=2, element_size=0.1)
       p0=Point(0,0,0,0.1)
       p1=Point(10,0,0,0.2)
       p2=Point(10,10,0,0.3)
       p3=Point(0,10,0,0.4)
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

       l6.setElementDistribution(5,1.2,True)
       l5.setElementDistribution(5,0.8,False)

       cl=CurveLoop(l0,l1,l2,l3)
       h=CurveLoop(l4,l5,l6)

       s=PlaneSurface(cl,holes=[h])

       d.addItems(PropertySet("test0", s))

       scrpt=d.getScriptString()
       ref = \
"""// generated by esys.pycad
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-02};
Point(2) = {1.00000000000000e+01, 0.00000000000000e+00, 0.00000000000000e+00 , 2.00000000000000e-02};
Point(3) = {1.00000000000000e+01, 1.00000000000000e+01, 0.00000000000000e+00 , 3.00000000000000e-02};
Point(4) = {0.00000000000000e+00, 1.00000000000000e+01, 0.00000000000000e+00 , 4.00000000000000e-02};
Point(5) = {5.00000000000000e+00, 5.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-04};
Point(6) = {7.00000000000000e+00, 5.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-04};
Point(7) = {5.00000000000000e+00, 7.00000000000000e+00, 0.00000000000000e+00 , 1.00000000000000e-04};
Line(8) = {1, 2};
Line(9) = {2, 3};
Line(10) = {3, 4};
Line(11) = {4, 1};
Line(12) = {5, 6};
Line(13) = {6, 7};
Transfinite Line{13} = 5 Using Progression 0.8;
Line(14) = {7, 5};
Transfinite Line{14} = 5 Using Bump 1.2;
Line Loop(15) = {8, 9, 10, 11};
Line Loop(16) = {12, 13, 14};
Plane Surface(17) = {15, 16};
Physical Surface(18) = {17};
"""
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

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
General.Terminal = 1;
General.ExpertMode = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.SubdivisionAlgorithm = 0;
Mesh.Smoothing = 1;
Mesh.RandomFactor = 1.00000000000000e-09;
Mesh.Algorithm = 1; // = MeshAdapt
Mesh.Algorithm3D = 4; // = Frontal
Point(1) = {-2.00000000000000e+00, -2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(2) = {2.00000000000000e+00, -2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(3) = {-2.00000000000000e+00, 2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(4) = {2.00000000000000e+00, 2.00000000000000e+00, -2.00000000000000e+00 , 1.00000000000000e-03};
Point(5) = {-2.00000000000000e+00, -2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
Point(6) = {2.00000000000000e+00, -2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
Point(7) = {-2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
Point(8) = {2.00000000000000e+00, 2.00000000000000e+00, 2.00000000000000e+00 , 1.00000000000000e-03};
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
Point(34) = {-1.00000000000000e+00 , -1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(35) = {1.00000000000000e+00 , -1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(36) = {-1.00000000000000e+00 , 1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(37) = {1.00000000000000e+00 , 1.00000000000000e+00, -1.00000000000000e+00 , 1.00000000000000e-03 };
Point(38) = {-1.00000000000000e+00 , -1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
Point(39) = {1.00000000000000e+00 , -1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
Point(40) = {-1.00000000000000e+00 , 1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
Point(41) = {1.00000000000000e+00 , 1.00000000000000e+00, 1.00000000000000e+00 , 1.00000000000000e-03 };
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
       self.assertEqual(scrpt.replace(' ',''), ref.replace(' ',''))

   def test_layer_cake(self):
       dom=layer_cake(GMSHDesign(dim=3, element_size=0.1),100.0,100.0,[10.,40.,80.,100.,150.])
       self.assertTrue(isinstance(dom,GMSHDesign),\
                         "LayerCake return is not a domain gmsh.Design.")


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

