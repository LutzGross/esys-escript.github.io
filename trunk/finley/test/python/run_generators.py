# $Id$


"""
checks the mesh generators against the reference meshes in test_meshes
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

import sys
import os
import unittest
from esys.escript import *
from esys.finley import Interval,Rectangle,Brick,JoinFaces, ReadGmsh, ReadMesh

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")
FINLEY_WORKDIR_PATH=FINLEY_WORKDIR

TEST_FILE_EXT=".test"
class Test_Generators(unittest.TestCase):

   def checker(self,dom,reference):
      dom_file=os.path.join(FINLEY_WORKDIR_PATH,TEST_FILE_EXT)
      dom.write(dom_file)
# Uncomment this section to dump the files for regression testing
#      if True:
#         dom.write(os.path.join(FINLEY_TEST_MESH_PATH,reference))
      dom_string=open(dom_file).read().splitlines()
      ref_string=open(os.path.join(FINLEY_TEST_MESH_PATH,reference)).read().splitlines()
      self.failUnlessEqual(len(dom_string),len(ref_string),"number of lines in mesh files does not match reference")
      for l in range(1,len(ref_string)):
         line=dom_string[l].strip()
         if os.name == "nt":
              line=line.replace("e+00","e+0").replace("e-00","e-0")
         self.failUnlessEqual(line,ref_string[l].strip(),"line %d=""%s"" in mesh file does not match line ""%s"" in reference %s" %(l,ref_string[l].strip(),line,os.path.join(FINLEY_TEST_MESH_PATH,reference)))

   def test_hex_1D_order1(self):
      file="hex_1D_order1.msh"
      my_dom=Interval(1,1)
      self.checker(my_dom,file)

   def test_hex_1D_order1_onFace(self):
      file="hex_1D_order1_onFace.msh"
      my_dom=Interval(1,1,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_1D_order2(self):
      file="hex_1D_order2.msh"
      my_dom=Interval(1,2)
      self.checker(my_dom,file)

   def test_hex_1D_order2_onFace(self):
      file="hex_1D_order2_onFace.msh"
      my_dom=Interval(1,2,useElementsOnFace=1)
      self.checker(my_dom,file)

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
      file="hex_contact_2D_order1.msh"
      ms1=Rectangle(1,1,1,l1=0.5,useElementsOnFace=False)
      ms2=Rectangle(1,1,1,l1=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2])
      self.checker(my_dom,file)

   def test_hex_contact_2D_order1_onFace(self):
      file="hex_contact_2D_order1_onFace.msh"
      ms1=Rectangle(1,1,1,l1=0.5,useElementsOnFace=True)
      ms2=Rectangle(1,1,1,l1=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2])
      self.checker(my_dom,file)

   def test_hex_contact_2D_order2(self):
      file="hex_contact_2D_order2.msh"
      ms1=Rectangle(1,1,2,l1=0.5,useElementsOnFace=False)
      ms2=Rectangle(1,1,2,l1=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2])
      self.checker(my_dom,file)

   def test_hex_contact_2D_order2_onFace(self):
      file="hex_contact_2D_order2_onFace.msh"
      ms1=Rectangle(1,1,2,l1=0.5,useElementsOnFace=True)
      ms2=Rectangle(1,1,2,l1=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2])
      self.checker(my_dom,file)

   def test_hex_contact_3D_order1(self):
      file="hex_contact_3D_order1.msh"
      ms1=Brick(1,1,1,1,l2=0.5,useElementsOnFace=False)
      ms2=Brick(1,1,1,1,l2=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2])
      self.checker(my_dom,file)

   def test_hex_contact_3D_order1_onFace(self):
      file="hex_contact_3D_order1_onFace.msh"
      ms1=Brick(1,1,1,1,l2=0.5,useElementsOnFace=True)
      ms2=Brick(1,1,1,1,l2=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2])
      self.checker(my_dom,file)

   def test_hex_contact_3D_order2(self):
      file="hex_contact_3D_order2.msh"
      ms1=Brick(1,1,1,2,l2=0.5,useElementsOnFace=False)
      ms2=Brick(1,1,1,2,l2=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2])
      self.checker(my_dom,file)

   def test_hex_contact_3D_order2_onFace(self):
      file="hex_contact_3D_order2_onFace.msh"
      ms1=Brick(1,1,1,2,l2=0.5,useElementsOnFace=True)
      ms2=Brick(1,1,1,2,l2=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2])
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
         self.failUnlessEqual(line,ref_string[l].strip(),"line %d=""%s"" in mesh file does not match line ""%s"" in reference %s" %(l,ref_string[l].strip(),line,reference_file))

   def test_Tri3(self):
       file="tri3_gmsh.msh"
       ref ="tri3.fly"
       test = os.path.join(FINLEY_WORKDIR,"tri3_test.fly")
       dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),2)
       dom.write(test)
       self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tri6(self):
       file="tri6_gmsh.msh"
       ref="tri6.fly"
       test = os.path.join(FINLEY_WORKDIR,"tri8_test.fly")
       dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),2)
       dom.write(test)
       self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tet4(self):
       file="tet4_gmsh.msh"
       ref="tet4.fly"
       test = os.path.join(FINLEY_WORKDIR,"tet4_test.fly")
       dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),3)
       dom.write(test)
       self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tet(self):
       file="tet10_gmsh.msh"
       ref="tet10.fly"
       test = os.path.join(FINLEY_WORKDIR,"tet10_test.fly")
       dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),3)
       dom.write(test)
       self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

class Test_Reader(unittest.TestCase):
   def test_ReadWriteTagNames(self):
       file="hex_2D_order2.msh"
       test = os.path.join(FINLEY_WORKDIR,"test.fly")
       dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,file),3)
       insertTagNames(dom,A=1,B=2)
       dom.write(test)
       dom2 = ReadMesh(test,3)
       t=getTagNames(dom)
       self.failUnless(len(t)==2)
       self.failUnless("A" in t)
       self.failUnless("B" in t)
       self.failUnless(dom2.getTag("A") == 1)
       self.failUnless(dom2.getTag("B") == 2)
       self.failUnless(dom2.isValidTagName("A"))
       self.failUnless(dom2.isValidTagName("B"))
       

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_Generators))
   suite.addTest(unittest.makeSuite(Test_GMSHReader))
   suite.addTest(unittest.makeSuite(Test_Reader))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
