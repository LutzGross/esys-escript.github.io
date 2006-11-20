# $Id: run_visualization_interface.py 798 2006-08-04 01:05:36Z gross $

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import sys
import unittest
import numarray
from esys.pycad import *

try:
     PYCAD_TEST_DATA=os.environ['PYCAD_TEST_DATA']
except KeyError:
     PYCAD_TEST_DATA='.'

try:
     PYCAD_WORKDIR=os.environ['PYCAD_WORKDIR']
except KeyError:
     PYCAD_WORKDIR='.'

PYCAD_TEST_MESH_PATH=PYCAD_TEST_DATA+"/data_meshes/"
PYCAD_WORKDIR_PATH=PYCAD_WORKDIR+"/"

class Test_PyCAD(unittest.TestCase):
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
if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_PyCAD))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
