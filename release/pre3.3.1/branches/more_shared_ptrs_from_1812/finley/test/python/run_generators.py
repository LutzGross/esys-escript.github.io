
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
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
checks the mesh generators against the reference meshes in test_meshes and test the finley integration schemes.
"""

import sys
import os
import unittest
from esys.escript import *
from esys.finley import Rectangle,Brick,JoinFaces, ReadGmsh, ReadMesh

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

FINLEY_TEST_MESH_PATH=FINLEY_TEST_DATA+"/data_meshes/"
#if os.name == "nt":
#   FINLEY_TEST_MESH_PATH = FINLEY_TEST_MESH_PATH+"win32/"
FINLEY_WORKDIR_PATH=FINLEY_WORKDIR+"/"

TEST_FILE_EXT=".test"
class Test_Generators(unittest.TestCase):

   def checker(self,dom,reference):
      if getMPISizeWorld() > 1: return
      dom_file=FINLEY_WORKDIR_PATH+TEST_FILE_EXT
      dom.write(dom_file)
# Uncomment this section to dump the files for regression testing
#      if True:
#         dom.write(FINLEY_TEST_MESH_PATH+reference)
      dom_string=open(dom_file).read().splitlines()
      ref_string=open(FINLEY_TEST_MESH_PATH+reference).read().splitlines()
      self.failUnlessEqual(len(dom_string),len(ref_string),"number of lines in mesh files does not match reference")
      for l in range(1,len(ref_string)):
	 line=dom_string[l].strip()
	 if os.name == "nt":
	       line=line.replace("e+00","e+0").replace("e-00","e-0")
         self.failUnlessEqual(line,ref_string[l].strip(),"line %d (%s) in mesh files does not match reference (%s)"%(l,ref_string[l].strip(),line))

   def test_hex_2D_order1(self):
      file="hex_2D_order1.msh"
      my_dom=Rectangle(1,1,1)
      self.checker(my_dom,file)

   def test_hex_2D_order1_onFace(self):
      file="hex_2D_order1_onFace.msh"
      my_dom=Rectangle(1,1,1,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_2D_order2(self):
      file="hex_2D_order2.msh"
      my_dom=Rectangle(1,1,2)
      self.checker(my_dom,file)

   def test_hex_2D_order2_onFace(self):
      file="hex_2D_order2_onFace.msh"
      my_dom=Rectangle(1,1,2,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_3D_order1(self):
      file="hex_3D_order1.msh"
      my_dom=Brick(1,1,1,1)
      self.checker(my_dom,file)

   def test_hex_3D_order1_onFace(self):
      file="hex_3D_order1_onFace.msh"
      my_dom=Brick(1,1,1,1,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_3D_order2(self):
      file="hex_3D_order2.msh"
      my_dom=Brick(1,1,1,2)
      self.checker(my_dom,file)

   def test_hex_3D_order2_onFace(self):
      file="hex_3D_order2_onFace.msh"
      my_dom=Brick(1,1,1,2,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_contact_2D_order1(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_2D_order1.msh"
      ms1=Rectangle(1,1,1,l1=0.5,useElementsOnFace=False)
      ms2=Rectangle(1,1,1,l1=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_2D_order1_onFace(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_2D_order1_onFace.msh"
      ms1=Rectangle(1,1,1,l1=0.5,useElementsOnFace=True)
      ms2=Rectangle(1,1,1,l1=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_2D_order2(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_2D_order2.msh"
      ms1=Rectangle(1,1,2,l1=0.5,useElementsOnFace=False)
      ms2=Rectangle(1,1,2,l1=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_2D_order2_onFace(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_2D_order2_onFace.msh"
      ms1=Rectangle(1,1,2,l1=0.5,useElementsOnFace=True)
      ms2=Rectangle(1,1,2,l1=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order1(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_3D_order1.msh"
      ms1=Brick(1,1,1,1,l2=0.5,useElementsOnFace=False)
      ms2=Brick(1,1,1,1,l2=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order1_onFace(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_3D_order1_onFace.msh"
      ms1=Brick(1,1,1,1,l2=0.5,useElementsOnFace=True)
      ms2=Brick(1,1,1,1,l2=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order2(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_3D_order2.msh"
      ms1=Brick(1,1,1,2,l2=0.5,useElementsOnFace=False)
      ms2=Brick(1,1,1,2,l2=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order2_onFace(self):
      if getMPISizeWorld() != 1: return
      file="hex_contact_3D_order2_onFace.msh"
      ms1=Brick(1,1,1,2,l2=0.5,useElementsOnFace=True)
      ms2=Brick(1,1,1,2,l2=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

class Test_GMSHReader(unittest.TestCase):
   def compare(self, test_file, reference_file):
      dom_string=open(test_file).read().splitlines()
      ref_string=open(reference_file).read().splitlines()
      self.failUnlessEqual(len(dom_string),len(ref_string),"number of lines in mesh files does not match reference")
      for l in range(1,len(ref_string)):
	 line=dom_string[l].strip()
	 if os.name == "nt":
	       line=line.replace("e+00","e+0").replace("e-00","e-0")
         self.failUnlessEqual(line,ref_string[l].strip(),"line %d (%s) in mesh files does not match reference (%s)"%(l,ref_string[l].strip(),line))

   def test_Tri3(self):
       # ReadGmsh is not MPI parallel
       if getMPISizeWorld() == 1:
         file="tri3_gmsh.msh"
         ref ="tri3.fly"
         test = FINLEY_WORKDIR+os.sep+"tri3_test.fly"
         dom = ReadGmsh(FINLEY_TEST_MESH_PATH+os.sep+file,2,optimize=False)
         dom.write(test)
         self.compare(test, FINLEY_TEST_MESH_PATH+os.sep+ref)

   def test_Tri6(self):
       # ReadGmsh is not MPI parallel
       if getMPISizeWorld() == 1:
         file="tri6_gmsh.msh"
         ref="tri6.fly"
         test = FINLEY_WORKDIR+os.sep+"tri8_test.fly"
         dom = ReadGmsh(FINLEY_TEST_MESH_PATH+os.sep+file,2,optimize=False)
         dom.write(test)
         self.compare(test, FINLEY_TEST_MESH_PATH+os.sep+ref)

   def test_Tet4(self):
       # ReadGmsh is not MPI parallel
       if getMPISizeWorld() == 1:
         file="tet4_gmsh.msh"
         ref="tet4.fly"
         test = FINLEY_WORKDIR+os.sep+"tet4_test.fly"
         dom = ReadGmsh(FINLEY_TEST_MESH_PATH+os.sep+file,3,optimize=False)
         dom.write(test)
         self.compare(test, FINLEY_TEST_MESH_PATH+os.sep+ref)

   def test_Tet10(self):
       # ReadGmsh is not MPI parallel
       if getMPISizeWorld() == 1:
         file="tet10_gmsh.msh"
         ref="tet10.fly"
         test = FINLEY_WORKDIR+os.sep+"tet10_test.fly"
         dom = ReadGmsh(FINLEY_TEST_MESH_PATH+os.sep+file,3,optimize=False)
         dom.write(test)
         self.compare(test, FINLEY_TEST_MESH_PATH+os.sep+ref)

class Test_Reader(unittest.TestCase):
   def test_ReadWriteTagNames(self):
       if getMPISizeWorld() != 1: return
       file="hex_2D_order2.msh"
       test = FINLEY_WORKDIR+os.sep+"test.fly"
       dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+file,3,optimize=False)
       insertTagNames(dom,A=1,B=2)
       dom.write(test)
       dom2 = ReadMesh(test,3,optimize=False)
       t=getTagNames(dom)
       self.failUnless(len(t)==6)
       self.failUnless("A" in t)
       self.failUnless("B" in t)
       self.failUnless(dom2.getTag("A") == 1)
       self.failUnless(dom2.getTag("B") == 2)
       self.failUnless(dom2.isValidTagName("A"))
       self.failUnless(dom2.isValidTagName("B"))

class Test_Integration(unittest.TestCase):
   TOL=EPSILON*100.
   def __test_2DQ(self,dom,order):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in xrange(order+1):
         for j in xrange(order+1):
             res=integrate(x[0]**i*x[1]**j)
             ref=1./((i+1)*(j+1))
             error=abs(res-ref)/abs(ref)
             self.failUnless(error<=self.TOL,"integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))

             res=integrate(x_bound[0]**i*x_bound[1]**j)
             ref=0
             if i == 0:
                ref += 2./(j+1)
             else:
                ref += 1./(j+1)
             if j == 0:
                ref += 2./(i+1)
             else:
                ref += 1./(i+1)
             error=abs(res-ref)/abs(ref)
             self.failUnless(error<=self.TOL,"surface integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))
            
   def __test_2DT(self,dom,order,raise_tol=1.):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in xrange(order+1):
         for j in xrange(order+1):
           if i+j<=order:
             res=integrate(x[0]**i*x[1]**j)
             ref=1./((i+1)*(j+1))
             error=abs(res-ref)/abs(ref)
             # print order,error
             self.failUnless(error<=self.TOL*raise_tol,"integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))

             res=integrate(x_bound[0]**i*x_bound[1]**j)
             ref=0
             if i == 0:
                ref += 2./(j+1)
             else:
                ref += 1./(j+1)
             if j == 0:
                ref += 2./(i+1)
             else:
                ref += 1./(i+1)
             error=abs(res-ref)/abs(ref)
             self.failUnless(error<=self.TOL*raise_tol,"surface integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))


   def __test_3DQ(self,dom,order):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in xrange(order+1):
         for j in xrange(order+1):
           for k in xrange(order+1):
             res=integrate(x[0]**i*x[1]**j*x[2]**k)
             ref=1./((i+1)*(j+1)*(k+1))
             error=abs(res-ref)/abs(ref)
             self.failUnless(error<=self.TOL,"integration for order (%s,%s,%s) failed. True value = %s, calculated = %s"%(i,j,k,ref,res))

             res=integrate(x_bound[0]**i*x_bound[1]**j*x_bound[2]**k)
             ref=0
             if i == 0:
                ref += 2./((j+1)*(k+1))
             else:
                ref += 1./((j+1)*(k+1))
             if j == 0:
                ref += 2./((i+1)*(k+1))
             else:
                ref += 1./((i+1)*(k+1))
             if k == 0:
                ref += 2./((i+1)*(j+1))
             else:
                ref += 1./((i+1)*(j+1))
             error=abs(res-ref)/abs(ref)
             self.failUnless(error<=self.TOL,"surface integration for order (%s,%s,%s) failed. True value = %s, calculated = %s"%(i,j,k,ref,res))

   def __test_3DT(self,dom,order,raise_tol=1.):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in xrange(order+1):
         for j in xrange(order+1):
           for k in xrange(order+1):
             if i+j+k<=order:
                res=integrate(x[0]**i*x[1]**j*x[2]**k)
                ref=1./((i+1)*(j+1)*(k+1))
                error=abs(res-ref)/abs(ref)
                self.failUnless(error<=self.TOL*raise_tol,"integration for order (%s,%s,%s) failed. True value = %s, calculated = %s"%(i,j,k,ref,res))

                res=integrate(x_bound[0]**i*x_bound[1]**j*x_bound[2]**k)
                ref=0
                if i == 0:
                   ref += 2./((j+1)*(k+1))
                else:
                   ref += 1./((j+1)*(k+1))
                if j == 0:
                   ref += 2./((i+1)*(k+1))
                else:
                   ref += 1./((i+1)*(k+1))
                if k == 0:
                   ref += 2./((i+1)*(j+1))
                else:
                   ref += 1./((i+1)*(j+1))
                error=abs(res-ref)/abs(ref)
                self.failUnless(error<=self.TOL*raise_tol,"surface integration for order (%s,%s,%s) failed. True value = %s, calculated = %s"%(i,j,k,ref,res))

   def test_hex2D_order1(self):
      my_dom=Rectangle(1,1,integrationOrder=1)
      self.__test_2DQ(my_dom,1)
   def test_hex2D_order2(self):
      my_dom=Rectangle(1,1,integrationOrder=2)
      self.__test_2DQ(my_dom,2)
   def test_hex2D_order3(self):
      my_dom=Rectangle(1,1,integrationOrder=3)
      self.__test_2DQ(my_dom,3)
   def test_hex2D_order4(self):
      my_dom=Rectangle(1,1,integrationOrder=4)
      self.__test_2DQ(my_dom,4)
   def test_hex2D_order5(self):
      my_dom=Rectangle(1,1,integrationOrder=5)
      self.__test_2DQ(my_dom,5)
   def test_hex2D_order6(self):
      my_dom=Rectangle(1,1,integrationOrder=6)
      self.__test_2DQ(my_dom,6)
   def test_hex2D_order7(self):
      my_dom=Rectangle(1,1,integrationOrder=7)
      self.__test_2DQ(my_dom,7)
   def test_hex2D_order8(self):
      my_dom=Rectangle(1,1,integrationOrder=8)
      self.__test_2DQ(my_dom,8)
   def test_hex2D_order9(self):
      my_dom=Rectangle(1,1,integrationOrder=9)
      self.__test_2DQ(my_dom,9)
   def test_hex2D_order10(self):
      my_dom=Rectangle(1,1,integrationOrder=10)
      self.__test_2DQ(my_dom,10)

   def test_Tet2D_order1(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=1)
      self.__test_2DT(my_dom,1)
   def test_Tet2D_order2(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=2)
      self.__test_2DT(my_dom,2)
   def test_Tet2D_order3(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=3)
      self.__test_2DT(my_dom,3)
   def test_Tet2D_order4(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=4)
      self.__test_2DT(my_dom,4)
   def test_Tet2D_order5(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=5)
      self.__test_2DT(my_dom,5)
   def test_Tet2D_order6(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=6)
      self.__test_2DT(my_dom,6)
   def test_Tet2D_order7(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=7)
      self.__test_2DT(my_dom,7)
   def test_Tet2D_order8(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=8)
      self.__test_2DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet2D_order9(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=9)
      self.__test_2DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet2D_order10(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tri3.fly",optimize=False,integrationOrder=10)
      self.__test_2DT(my_dom,10)

   def test_hex3D_order1(self):
      my_dom=Brick(1,1,1,integrationOrder=1)
      self.__test_3DQ(my_dom,1)
   def test_hex3D_order2(self):
      my_dom=Brick(1,1,1,integrationOrder=2)
      self.__test_3DQ(my_dom,2)
   def test_hex3D_order3(self):
      my_dom=Brick(1,1,1,integrationOrder=3)
      self.__test_3DQ(my_dom,3)
   def test_hex3D_order4(self):
      my_dom=Brick(1,1,1,integrationOrder=4)
      self.__test_3DQ(my_dom,4)
   def test_hex3D_order5(self):
      my_dom=Brick(1,1,1,integrationOrder=5)
      self.__test_3DQ(my_dom,5)
   def test_hex3D_order6(self):
      my_dom=Brick(1,1,1,integrationOrder=6)
      self.__test_3DQ(my_dom,6)
   def test_hex3D_order7(self):
      my_dom=Brick(1,1,1,integrationOrder=7)
      self.__test_3DQ(my_dom,7)
   def test_hex3D_order8(self):
      my_dom=Brick(1,1,1,integrationOrder=8)
      self.__test_3DQ(my_dom,8)
   def test_hex3D_order9(self):
      my_dom=Brick(1,1,1,integrationOrder=9)
      self.__test_3DQ(my_dom,9)
   def test_hex3D_order10(self):
      my_dom=Brick(1,1,1,integrationOrder=10)
      self.__test_3DQ(my_dom,10)

   def test_Tet3D_order1(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=1)
      self.__test_3DT(my_dom,1)
   def test_Tet3D_order2(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=2)
      self.__test_3DT(my_dom,2)
   def test_Tet3D_order3(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=3)
      self.__test_3DT(my_dom,3)
   def test_Tet3D_order4(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=4)
      self.__test_3DT(my_dom,4)
   def test_Tet3D_order5(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=5)
      self.__test_3DT(my_dom,5)
   def test_Tet3D_order6(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=6)
      self.__test_3DT(my_dom,6,1./sqrt(EPSILON))
   def test_Tet3D_order7(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=7)
      self.__test_3DT(my_dom,7,1./sqrt(EPSILON))
   def test_Tet3D_order8(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=8)
      self.__test_3DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet3D_order9(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=9)
      self.__test_3DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet3D_order10(self):
      my_dom = ReadMesh(FINLEY_TEST_MESH_PATH+os.sep+"tet4.fly",optimize=False,integrationOrder=10)
      self.__test_3DT(my_dom,10)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_Generators))
   suite.addTest(unittest.makeSuite(Test_GMSHReader))
   suite.addTest(unittest.makeSuite(Test_Reader))
   suite.addTest(unittest.makeSuite(Test_Integration))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

