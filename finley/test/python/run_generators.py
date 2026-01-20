
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
checks the mesh generators against the reference meshes in test_meshes and test the finley integration schemes.
"""

import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.finley import Rectangle, Brick, JoinFaces, ReadGmsh, ReadMesh

mpisize = getMPISizeWorld()

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")
#if os.name == "nt":
#   FINLEY_TEST_MESH_PATH = os.path.join(FINLEY_TEST_MESH_PATH,"win32")
FINLEY_WORKDIR_PATH=FINLEY_WORKDIR

TEST_FILE_PRE="test_"

@unittest.skipIf(mpisize>1, "multiple processes not supported for mesh writes")
class Test_GeneratorsOnFinley(unittest.TestCase):

   def checker(self, dom, reference):
      dom_file=os.path.join(FINLEY_WORKDIR_PATH, TEST_FILE_PRE+reference)
      dom.write(dom_file)
# Uncomment this section to dump the files for regression testing
      #if True:
      #   dom.write(os.path.join(FINLEY_TEST_MESH_PATH,reference))
      dom_string=open(dom_file).read().splitlines()
      ref_string=open(os.path.join(FINLEY_TEST_MESH_PATH,reference)).read().splitlines()
      self.assertEqual(len(dom_string),len(ref_string),"number of lines in mesh files does not match reference")
      taglist=[]
      taglist_ref=[]
      reading_tags=False
      # compare files but ignore tags for now since they might be in different
      # order
      for l in range(1, len(ref_string)):
         line=dom_string[l].strip()
         if not reading_tags:
             self.assertEqual(line,ref_string[l].strip(),"line %d (%s) in mesh file does not match reference (%s)"%(l,ref_string[l].strip(),line))
         else:
             taglist.append(line)
             taglist_ref.append(ref_string[l].strip())
         if line=='Tags':
             reading_tags=True
      # now compare tag lists disregarding order
      self.assertEqual(sorted(taglist),sorted(taglist_ref),"list of tags in mesh file does not match reference.")

   def test_hex_2D_order1(self):
      file="hex_2D_order1.msh"
      my_dom=Rectangle(1,1,1,useElementsOnFace=0)
      self.checker(my_dom,file)

   def test_hex_2D_order1_onFace(self):
      file="hex_2D_order1_onFace.msh"
      my_dom=Rectangle(1,1,1,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_2D_order2(self):
      file="hex_2D_order2.msh"
      my_dom=Rectangle(1,1,2, useElementsOnFace=0)
      self.checker(my_dom,file)

   def test_hex_2D_order1_macro(self):
      file="hex_2D_order1_macro.msh"
      my_dom=Rectangle(1,1,-1,useElementsOnFace=0)
      self.checker(my_dom,file)

   def test_hex_2D_order2_onFace(self):
      file="hex_2D_order2_onFace.msh"
      my_dom=Rectangle(1,1,2,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_3D_order1(self):
      file="hex_3D_order1.msh"
      my_dom=Brick(1,1,1,1, useElementsOnFace=0)
      self.checker(my_dom,file)

   def test_hex_3D_order1_onFace(self):
      file="hex_3D_order1_onFace.msh"
      my_dom=Brick(1,1,1,1,useElementsOnFace=1)
      self.checker(my_dom,file)

   def test_hex_3D_order2(self):
      file="hex_3D_order2.msh"
      my_dom=Brick(1,1,1,2)
      self.checker(my_dom,file)

   def test_hex_3D_order2(self):
      file="hex_3D_order1_macro.msh"
      my_dom=Brick(1,1,1,-1, useElementsOnFace=0)
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
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_2D_order1_onFace(self):
      file="hex_contact_2D_order1_onFace.msh"
      ms1=Rectangle(1,1,1,l1=0.5,useElementsOnFace=True)
      ms2=Rectangle(1,1,1,l1=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_2D_order2(self):
      file="hex_contact_2D_order2.msh"
      ms1=Rectangle(1,1,2,l1=0.5,useElementsOnFace=False)
      ms2=Rectangle(1,1,2,l1=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_2D_order2_onFace(self):
      file="hex_contact_2D_order2_onFace.msh"
      ms1=Rectangle(1,1,2,l1=0.5,useElementsOnFace=True)
      ms2=Rectangle(1,1,2,l1=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order1(self):
      file="hex_contact_3D_order1.msh"
      ms1=Brick(1,1,1,1,l2=0.5,useElementsOnFace=False)
      ms2=Brick(1,1,1,1,l2=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order1_onFace(self):
      file="hex_contact_3D_order1_onFace.msh"
      ms1=Brick(1,1,1,1,l2=0.5,useElementsOnFace=True)
      ms2=Brick(1,1,1,1,l2=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order2(self):
      file="hex_contact_3D_order2.msh"
      ms1=Brick(1,1,1,2,l2=0.5,useElementsOnFace=False)
      ms2=Brick(1,1,1,2,l2=0.5,useElementsOnFace=False)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

   def test_hex_contact_3D_order2_onFace(self):
      file="hex_contact_3D_order2_onFace.msh"
      ms1=Brick(1,1,1,2,l2=0.5,useElementsOnFace=True)
      ms2=Brick(1,1,1,2,l2=0.5,useElementsOnFace=True)
      ms2.setX(ms2.getX()+[0,0,0.5])
      my_dom=JoinFaces([ms1,ms2],optimize=False)
      self.checker(my_dom,file)

@unittest.skipIf(mpisize>1, "multiple processes not supported for mesh writes")
class Test_GmshReaderOnFinley(unittest.TestCase):
   def compare(self, test_file, reference_file):
      dom_string=open(test_file).read().splitlines()
      ref_string=open(reference_file).read().splitlines()
      self.assertEqual(len(dom_string),len(ref_string),"number of lines in mesh files does not match reference")
      for l in range(1,len(ref_string)):
         line=dom_string[l].strip()
         self.assertEqual(line,ref_string[l].strip(),"line %d (%s) in mesh files does not match reference (%s)"%(l,ref_string[l].strip(),line))

   def test_Tri3(self):
         file="tri3_gmsh.msh"
         ref ="tri3.fly"
         test = os.path.join(FINLEY_WORKDIR,"tri3_test.fly")
         dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),2,optimize=False)
         dom.write(test)
         self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tri6(self):
         file="tri6_gmsh.msh"
         ref="tri6.fly"
         test = os.path.join(FINLEY_WORKDIR,"tri6_test.fly")
         dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),2,optimize=False)
         dom.write(test)
         self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tri6_macro(self):
         file="tri6_gmsh.msh"
         ref="tri6_macro.fly"
         test = os.path.join(FINLEY_WORKDIR,"tri6_macro_test.fly")
         dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),2,useMacroElements=True,optimize=False)
         dom.write(test)
         self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tet4(self):
         file="tet4_gmsh.msh"
         ref="tet4.fly"
         test = os.path.join(FINLEY_WORKDIR,"tet4_test.fly")
         dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),3,optimize=False)
         dom.write(test)
         self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tet10(self):
         file="tet10_gmsh.msh"
         ref="tet10.fly"
         test = os.path.join(FINLEY_WORKDIR,"tet10_test.fly")
         dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),3,optimize=False)
         dom.write(test)
         self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

   def test_Tet10_macro(self):
         file="tet10_gmsh.msh"
         ref="tet10_macro.fly"
         test = os.path.join(FINLEY_WORKDIR,"tet10_macro_test.fly")
         dom = ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH,file),3,useMacroElements=True,optimize=False)
         dom.write(test)
         self.compare(test, os.path.join(FINLEY_TEST_MESH_PATH,ref))

@unittest.skipIf(mpisize>1, "multiple processes not supported for mesh writes")
class Test_ReaderOnFinley(unittest.TestCase):
   def test_ReadWriteTagNames(self):
       file="hex_2D_order2.msh"
       test = os.path.join(FINLEY_WORKDIR,"test.fly")
       dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,file),3,optimize=False)
       insertTagNames(dom,A=1,B=2)
       dom.write(test)
       dom2 = ReadMesh(test,3,optimize=False)
       t=getTagNames(dom)
       self.assertTrue(len(t)==6)
       self.assertTrue("A" in t)
       self.assertTrue("B" in t)
       self.assertTrue(dom2.getTag("A") == 1)
       self.assertTrue(dom2.getTag("B") == 2)
       self.assertTrue(dom2.isValidTagName("A"))
       self.assertTrue(dom2.isValidTagName("B"))

class Test_IntegrationOnFinley(unittest.TestCase):
   TOL=EPSILON*500.
   def __test_2DQ(self,dom,order):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in range(order+1):
         for j in range(order+1):
             res=integrate(x[0]**i*x[1]**j)
             ref=1./((i+1)*(j+1))
             error=abs(res-ref)/abs(ref)
             self.assertTrue(error<=self.TOL,"integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))

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
             self.assertTrue(error<=self.TOL,"surface integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))

   def __test_2DT(self,dom,order,raise_tol=1.):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in range(order+1):
         for j in range(order+1):
           if i+j<=order:
             res=integrate(x[0]**i*x[1]**j)
             ref=1./((i+1)*(j+1))
             error=abs(res-ref)/abs(ref)
             # print order,error
             self.assertTrue(error<=self.TOL*raise_tol,"integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))

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
             self.assertTrue(error<=self.TOL*raise_tol,"surface integration for order (%s,%s) failed. True value = %s, calculated = %s"%(i,j,ref,res))


   def __test_3DQ(self,dom,order):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in range(order+1):
         for j in range(order+1):
           for k in range(order+1):
             res=integrate(x[0]**i*x[1]**j*x[2]**k)
             ref=1./((i+1)*(j+1)*(k+1))
             error=abs(res-ref)/abs(ref)
             self.assertTrue(error<=self.TOL,"integration for order (%s,%s,%s) failed. True value = %s, calculated = %s (error=%e)"%(i,j,k,ref,res, error))

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
             self.assertTrue(error<=self.TOL,"surface integration for order (%s,%s,%s) failed. True value = %s, calculated = %s (error=%e)"%(i,j,k,ref,res, error))

   def __test_3DT(self,dom,order,raise_tol=1.):
       x=Function(dom).getX()
       x_bound=FunctionOnBoundary(dom).getX()
       for i in range(order+1):
         for j in range(order+1):
           for k in range(order+1):
             if i+j+k<=order:
                res=integrate(x[0]**i*x[1]**j*x[2]**k)
                ref=1./((i+1)*(j+1)*(k+1))
                error=abs(res-ref)/abs(ref)
                self.assertTrue(error<=self.TOL*raise_tol,"integration for order (%s,%s,%s) failed. True value = %s, calculated = %s (error=%e)"%(i,j,k,ref,res,error))

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
                self.assertTrue(error<=self.TOL*raise_tol,"surface integration for order (%s,%s,%s) failed. True value = %s, calculated = %s (error=%e)"%(i,j,k,ref,res,error))

   #==========================================================================
   def test_hex2D_order1_integorder1(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=1)
      self.__test_2DQ(my_dom,1)
   def test_hex2D_order1_integorder2(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=2)
      self.__test_2DQ(my_dom,2)
   def test_hex2D_order1_integorder3(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=3)
      self.__test_2DQ(my_dom,3)
   def test_hex2D_order1_integorder4(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=4)
      self.__test_2DQ(my_dom,4)
   def test_hex2D_order1_integorder5(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=5)
      self.__test_2DQ(my_dom,5)
   def test_hex2D_order1_integorder6(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=6)
      self.__test_2DQ(my_dom,6)
   def test_hex2D_order1_integorder7(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=7)
      self.__test_2DQ(my_dom,7)
   def test_hex2D_order1_integorder8(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=8)
      self.__test_2DQ(my_dom,8)
   def test_hex2D_order1_integorder9(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=9)
      self.__test_2DQ(my_dom,9)
   def test_hex2D_order1_integorder10(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,integrationOrder=10)
      self.__test_2DQ(my_dom,10)
   #==========================================================================
   def test_hex2D_order2_integorder1(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=1)
      self.__test_2DQ(my_dom,1)
   def test_hex2D_order2_integorder2(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=2)
      self.__test_2DQ(my_dom,2)
   def test_hex2D_order2_integorder3(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=3)
      self.__test_2DQ(my_dom,3)
   def test_hex2D_order2_integorder4(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=4)
      self.__test_2DQ(my_dom,4)
   def test_hex2D_order2_integorder5(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=5)
      self.__test_2DQ(my_dom,5)
   def test_hex2D_order2_integorder6(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=6)
      self.__test_2DQ(my_dom,6)
   def test_hex2D_order2_integorder7(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=7)
      self.__test_2DQ(my_dom,7)
   def test_hex2D_order2_integorder8(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=8)
      self.__test_2DQ(my_dom,8)
   def test_hex2D_order2_integorder9(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=9)
      self.__test_2DQ(my_dom,9)
   def test_hex2D_order2_integorder10(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=2,integrationOrder=10)
      self.__test_2DQ(my_dom,10)
   #==========================================================================
   def test_hex2D_macro_integorder1(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=1)
      self.__test_2DQ(my_dom,1)
   def test_hex2D_macro_integmacro(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=2)
      self.__test_2DQ(my_dom,2)
   def test_hex2D_macro_integorder3(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=3)
      self.__test_2DQ(my_dom,3)
   def test_hex2D_macro_integorder4(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=4)
      self.__test_2DQ(my_dom,4)
   def test_hex2D_macro_integorder5(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=5)
      self.__test_2DQ(my_dom,5)
   def test_hex2D_macro_integorder6(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=6)
      self.__test_2DQ(my_dom,6)
   def test_hex2D_macro_integorder7(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=7)
      self.__test_2DQ(my_dom,7)
   def test_hex2D_macro_integorder8(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=8)
      self.__test_2DQ(my_dom,8)
   def test_hex2D_macro_integorder9(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=9)
      self.__test_2DQ(my_dom,9)
   def test_hex2D_macro_integorder10(self):
      NE=getMPIRankWorld()
      my_dom=Rectangle(NE,NE,order=-1,useElementsOnFace=0,integrationOrder=10)
      self.__test_2DQ(my_dom,10)
   #==========================================================================
   def test_Tet2D_order1_integorder1(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=1)
      self.__test_2DT(my_dom,1)
   def test_Tet2D_order1_integorder2(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=2)
      self.__test_2DT(my_dom,2)
   def test_Tet2D_order1_integorder3(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=3)
      self.__test_2DT(my_dom,3)
   def test_Tet2D_order1_integorder4(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=4)
      self.__test_2DT(my_dom,4)
   def test_Tet2D_order1_integorder5(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=5)
      self.__test_2DT(my_dom,5)
   def test_Tet2D_order1_integorder6(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=6)
      self.__test_2DT(my_dom,6)
   def test_Tet2D_order1_integorder7(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=7)
      self.__test_2DT(my_dom,7)
   def test_Tet2D_order1_integorder8(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=8)
      self.__test_2DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet2D_order1_integorder9(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=9)
      self.__test_2DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet2D_order1_integorder10(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri3.fly"),optimize=False,integrationOrder=10)
      self.__test_2DT(my_dom,10)
   #==========================================================================
   def test_Tet2D_order2_integorder2(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=1)
      self.__test_2DT(my_dom,1)
   def test_Tet2D_order2_integorder2(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=2)
      self.__test_2DT(my_dom,2)
   def test_Tet2D_order2_integorder3(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=3)
      self.__test_2DT(my_dom,3)
   def test_Tet2D_order2_integorder4(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=4)
      self.__test_2DT(my_dom,4)
   def test_Tet2D_order2_integorder5(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=5)
      self.__test_2DT(my_dom,5)
   def test_Tet2D_order2_integorder6(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=6)
      self.__test_2DT(my_dom,6)
   def test_Tet2D_order2_integorder7(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=7)
      self.__test_2DT(my_dom,7)
   def test_Tet2D_order2_integorder8(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=8)
      self.__test_2DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet2D_order2_integorder9(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=9)
      self.__test_2DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet2D_order2_integorder10(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6.fly"),optimize=False,integrationOrder=10)
      self.__test_2DT(my_dom,10)
   #==========================================================================
   def test_Tet2D_macro_integmacro(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=1)
      self.__test_2DT(my_dom,1)
   def test_Tet2D_macro_integmacro(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=2)
      self.__test_2DT(my_dom,2)
   def test_Tet2D_macro_integorder3(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=3)
      self.__test_2DT(my_dom,3)
   def test_Tet2D_macro_integorder4(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=4)
      self.__test_2DT(my_dom,4)
   def test_Tet2D_macro_integorder5(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=5)
      self.__test_2DT(my_dom,5)
   def test_Tet2D_macro_integorder6(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=6)
      self.__test_2DT(my_dom,6)
   def test_Tet2D_macro_integorder7(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=7)
      self.__test_2DT(my_dom,7)
   def test_Tet2D_macro_integorder8(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=8)
      self.__test_2DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet2D_macro_integorder9(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=9)
      self.__test_2DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet2D_macro_integorder10(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tri6_macro.fly"),optimize=False,integrationOrder=10)
      self.__test_2DT(my_dom,10)
   #==========================================================================
   def test_hex3D_order1_integorder1(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=1)
      self.__test_3DQ(my_dom,1)
   def test_hex3D_order1_integorder2(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=2)
      self.__test_3DQ(my_dom,2)
   def test_hex3D_order1_integorder3(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=3)
      self.__test_3DQ(my_dom,3)
   def test_hex3D_order1_integorder4(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=4)
      self.__test_3DQ(my_dom,4)
   def test_hex3D_order1_integorder5(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=5)
      self.__test_3DQ(my_dom,5)
   def test_hex3D_order1_integorder6(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=6)
      self.__test_3DQ(my_dom,6)
   def test_hex3D_order1_integorder7(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=7)
      self.__test_3DQ(my_dom,7)
   def test_hex3D_order1_integorder8(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=8)
      self.__test_3DQ(my_dom,8)
   def test_hex3D_order1_integorder9(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=9)
      self.__test_3DQ(my_dom,9)
   def test_hex3D_order1_integorder10(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,integrationOrder=10)
      self.__test_3DQ(my_dom,10)
   #==========================================================================
   def test_hex3D_order2_integorder2(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=1)
      self.__test_3DQ(my_dom,1)
   def test_hex3D_order2_integorder2(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=2)
      self.__test_3DQ(my_dom,2)
   def test_hex3D_order2_integorder3(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=3)
      self.__test_3DQ(my_dom,3)
   def test_hex3D_order2_integorder4(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=4)
      self.__test_3DQ(my_dom,4)
   def test_hex3D_order2_integorder5(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=5)
      self.__test_3DQ(my_dom,5)
   def test_hex3D_order2_integorder6(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=6)
      self.__test_3DQ(my_dom,6)
   def test_hex3D_order2_integorder7(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=7)
      self.__test_3DQ(my_dom,7)
   def test_hex3D_order2_integorder8(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=8)
      self.__test_3DQ(my_dom,8)
   def test_hex3D_order2_integorder9(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=9)
      self.__test_3DQ(my_dom,9)
   def test_hex3D_order2_integorder10(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=2,integrationOrder=10)
      self.__test_3DQ(my_dom,10)
   #==========================================================================
   def test_hex3D_macro_integmacro(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,useElementsOnFace=0,order=-1,integrationOrder=1)
      self.__test_3DQ(my_dom,1)
   def test_hex3D_macro_integmacro(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=2)
      self.__test_3DQ(my_dom,2)
   def test_hex3D_macro_integorder3(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=3)
      self.__test_3DQ(my_dom,3)
   def test_hex3D_macro_integorder4(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=4)
      self.__test_3DQ(my_dom,4)
   def test_hex3D_macro_integorder5(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=5)
      self.__test_3DQ(my_dom,5)
   def test_hex3D_macro_integorder6(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=6)
      self.__test_3DQ(my_dom,6)
   def test_hex3D_macro_integorder7(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=7)
      self.__test_3DQ(my_dom,7)
   def test_hex3D_macro_integorder8(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=8)
      self.__test_3DQ(my_dom,8)
   def test_hex3D_macro_integorder9(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1, useElementsOnFace=0,integrationOrder=9)
      self.__test_3DQ(my_dom,9)
   def test_hex3D_macro_integorder10(self):
      NE=getMPIRankWorld()
      my_dom=Brick(NE,NE,NE,order=-1,useElementsOnFace=0,integrationOrder=10)
      self.__test_3DQ(my_dom,10)
   #==========================================================================
   def test_Tet3D_order1_integorder1(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=1)
      self.__test_3DT(my_dom,1)
   def test_Tet3D_order1_integorder2(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=2)
      self.__test_3DT(my_dom,2)
   def test_Tet3D_order1_integorder3(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=3)
      self.__test_3DT(my_dom,3)
   def test_Tet3D_order1_integorder4(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=4)
      self.__test_3DT(my_dom,4)
   def test_Tet3D_order1_integorder5(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=5)
      self.__test_3DT(my_dom,5)
   def test_Tet3D_order1_integorder6(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=6)
      self.__test_3DT(my_dom,6,1./sqrt(EPSILON))
   def test_Tet3D_order1_integorder7(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=7)
      self.__test_3DT(my_dom,7,1./sqrt(EPSILON))
   def test_Tet3D_order1_integorder8(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=8)
      self.__test_3DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet3D_order1_integorder9(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=9)
      self.__test_3DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet3D_order1_integorder10(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet4.fly"),optimize=False,integrationOrder=10)
      self.__test_3DT(my_dom,10,1./sqrt(EPSILON))
   #==========================================================================
   def test_Tet3D_order2_integorder2(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=1)
      self.__test_3DT(my_dom,1)
   def test_Tet3D_order2_integorder2(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=2)
      self.__test_3DT(my_dom,2)
   def test_Tet3D_order2_integorder3(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=3)
      self.__test_3DT(my_dom,3)
   def test_Tet3D_order2_integorder4(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=4)
      self.__test_3DT(my_dom,4)
   def test_Tet3D_order2_integorder5(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=5)
      self.__test_3DT(my_dom,5)
   def test_Tet3D_order2_integorder6(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=6)
      self.__test_3DT(my_dom,6,1./sqrt(EPSILON))
   def test_Tet3D_order2_integorder7(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=7)
      self.__test_3DT(my_dom,7,1./sqrt(EPSILON))
   def test_Tet3D_order2_integorder8(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=8)
      self.__test_3DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet3D_order2_integorder9(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=9)
      self.__test_3DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet3D_order2_integorder10(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10.fly"),optimize=False,integrationOrder=10)
      self.__test_3DT(my_dom,10,1./sqrt(EPSILON))
   #==========================================================================
   def test_Tet3D_macro_integmacro(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=1)
      self.__test_3DT(my_dom,1)
   def test_Tet3D_macro_integmacro(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=2)
      self.__test_3DT(my_dom,2)
   def test_Tet3D_macro_integorder3(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=3)
      self.__test_3DT(my_dom,3)
   def test_Tet3D_macro_integorder4(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=4)
      self.__test_3DT(my_dom,4)
   def test_Tet3D_macro_integorder5(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=5)
      self.__test_3DT(my_dom,5)
   def test_Tet3D_macro_integorder6(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=6)
      self.__test_3DT(my_dom,6,1./sqrt(EPSILON))
   def test_Tet3D_macro_integorder7(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=7)
      self.__test_3DT(my_dom,7,1./sqrt(EPSILON))
   def test_Tet3D_macro_integorder8(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=8)
      self.__test_3DT(my_dom,8,1./sqrt(EPSILON))
   def test_Tet3D_macro_integorder9(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=9)
      self.__test_3DT(my_dom,9,1./sqrt(EPSILON))
   def test_Tet3D_macro_integorder10(self):
      my_dom = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"tet10_macro.fly"),optimize=False,integrationOrder=10)
      self.__test_3DT(my_dom,10,1./sqrt(EPSILON))

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

