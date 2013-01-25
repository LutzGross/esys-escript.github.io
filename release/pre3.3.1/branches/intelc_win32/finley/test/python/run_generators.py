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
from esys.finley import Interval,Rectangle,Brick,JoinFaces

FINLEY_TEST_DATA_PATH=os.environ['FINLEY_TEST_DATA']+"/data_meshes/"
if os.name == "nt":
   FINLEY_TEST_DATA_PATH = FINLEY_TEST_DATA_PATH+"win32/"
FINLEY_WORKDIR_PATH=os.environ['FINLEY_WORKDIR']+"/"

TEST_FILE_EXT=".test"
class Test_Generators(unittest.TestCase):

   def checker(self,dom,reference):
      dom_file=FINLEY_WORKDIR_PATH+TEST_FILE_EXT
      dom.write(dom_file)
# Uncomment this section to dump the files for regression testing
#      if True:
#         dom.write(FINLEY_TEST_DATA_PATH+reference)
      dom_string=open(dom_file).read().splitlines()
      ref_string=open(FINLEY_TEST_DATA_PATH+reference).read().splitlines()
      self.failUnlessEqual(len(dom_string),len(ref_string),"number of lines in mesh files does not match reference")
      for l in range(1,len(ref_string)):
         self.failUnlessEqual(dom_string[l].strip(),ref_string[l].strip(),"line %d (%s) in mesh files does not match reference (%s)"%(l,ref_string[l].strip(),dom_string[l].strip()))

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

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_Generators))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
