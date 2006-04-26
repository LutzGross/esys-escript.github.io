# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import sys
import unittest
from esys.escript import *
from esys.finley import ReadMesh

FINLEY_TEST_MESH_PATH=os.environ['FINLEY_TEST_DATA']+"/data_meshes/"
FINLEY_WORKDIR_PATH=os.environ['FINLEY_WORKDIR']+"/"

class Test_VisualizationInterface(unittest.TestCase):
   def check_vtk(self,f,reference_f):
      out_string=open(FINLEY_WORKDIR_PATH+f).read().splitlines()
      ref_string=open(FINLEY_TEST_MESH_PATH+reference_f).read().splitlines()
      c=0
      for l in range(0,len(ref_string)):
         if not ref_string[l].strip()[:2]=="<!":
           self.failUnlessEqual(out_string[c].strip().replace("e-00","e+00"),ref_string[l].strip(),"line %d (%s) in vtk files does not match reference (%s)"%(c,out_string[c].strip(),ref_string[l].strip()))
           c+=1

   def check_dx(self,f,reference_f):
      out_string=open(FINLEY_WORKDIR_PATH+f).read().splitlines()
      ref_string=open(FINLEY_TEST_MESH_PATH+reference_f).read().splitlines()
      c=0
      for l in range(0,len(ref_string)):
         if not ref_string[l].strip()[0]=="#":
           self.failUnlessEqual(out_string[c].strip().replace("e-00","e+00"),ref_string[l].strip(),"line %d (%s) in dx file does not match reference (%s)"%(c,out_string[c].strip(),ref_string[l].strip()))
           c+=1

class Test_VTKFiles(Test_VisualizationInterface):
  # ======================================================================================================================
  def test_hex_2D_order2_vtk(self):
     reference="hex_2D_o2.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     saveVTK(FINLEY_WORKDIR_PATH+"hex_2D_order2.xml",domain=dom)
     self.check_vtk("hex_2D_order2.xml",reference)

  def test_hex_2D_order2_AllPoints_Scalar_vtk(self):
     reference="hex_2D_o1_node_3xs.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_2D_order2_AllPoints_Scalar.xml",data_r=x_r[0],data_n=x_n[0],data=x[0])
     self.check_vtk("hex_2D_order2_AllPoints_Scalar.xml",reference)
  def test_hex_2D_order2_02Points_Scalar_vtk(self):
     reference="hex_2D_o2_node_2xs.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Solution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_2D_order2_O2Points_Scalar.xml",data_n=x_n[0],data=x[0])
     self.check_vtk("hex_2D_order2_O2Points_Scalar.xml",reference)
  def test_hex_2D_order2_2Cells_Scalar_vtk(self):
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     try: 
        saveVTK(FINLEY_WORKDIR_PATH+"hex_2D_order2_2Cells_Scalar.xml",data=x[0],data_b=x_b[0])
        self.fail("non-matching data not detected.")
     except StandardError:
        pass
  def test_hex_2D_order2_BoundrayPoint_Scalar_vtk(self):
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     try: 
        saveVTK(FINLEY_WORKDIR_PATH+"hex_2D_order2_BoundrayPoint_Scalar.xml",data=x[0],data_b=x_b[0])
        self.fail("non-matching data not detected.")
     except StandardError:
        pass
  def test_hex_2D_order2_Cells_AllData_vtk(self):
     reference="hex_2D_o2_cell_all.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_2D_order2_Cells_AllData.xml",data_s=x[0],data_v=x[0]*[1.,2.],data_t=x[0]*[[11.,12.],[21.,22.]],data_t2=x[0]*[[-11.,-12.],[-21.,-22.]])
     self.check_vtk("hex_2D_order2_Cells_AllData.xml",reference)

  def test_hex_2D_order2_CellsPoints_AllData_vtk(self):
     reference="hex_2D_o2_cellnode_all.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_2D_order2_CellsPoints_AllData.xml",data_sp=x_p[0],
                                                     data_vp=x_p[0]*[1.,2.],
                                                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                                                     data_sc=x_c[0],
                                                     data_vc=x_c[0]*[1.,2.],
                                                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_2D_order2_CellsPoints_AllData.xml",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order1_ContinuousFunction_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_2D_order1_Solution_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Solution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_Solution_Scalar.xml",reference)
  def test_hex_contact_2D_order1_Solution_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Solution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_Solution_Vector.xml",reference)
  def test_hex_contact_2D_order1_Solution_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Solution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_Solution_Tensor.xml",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ReducedSolution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_2D_order1_Function_Scalar_vtk(self):
     reference="hex_2D_o1_cell_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Function_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_Function_Scalar.xml",reference)
  def test_hex_contact_2D_order1_Function_Vector_vtk(self):
     reference="hex_2D_o1_cell_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Function_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_Function_Vector.xml",reference)
  def test_hex_contact_2D_order1_Function_Tensor_vtk(self):
     reference="hex_2D_o1_cell_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Function_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_Function_Tensor.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o1_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o1_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o1_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o1_f_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o1_f_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o1_f_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order1_onFace_FunctionOnContactOne_Tensor.xml",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order2_ContinuousFunction_Scalar_vtk(self):
     reference="hex_2D_o2_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Vector_vtk(self):
     reference="hex_2D_o2_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Tensor_vtk(self):
     reference="hex_2D_o2_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_2D_order2_Solution_Scalar_vtk(self):
     reference="hex_2D_o2_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Solution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_Solution_Scalar.xml",reference)
  def test_hex_contact_2D_order2_Solution_Vector_vtk(self):
     reference="hex_2D_o2_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Solution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_Solution_Vector.xml",reference)
  def test_hex_contact_2D_order2_Solution_Tensor_vtk(self):
     reference="hex_2D_o2_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Solution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_Solution_Tensor.xml",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ReducedSolution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_2D_order2_Function_Scalar_vtk(self):
     reference="hex_2D_o2_cell_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Function_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_Function_Scalar.xml",reference)
  def test_hex_contact_2D_order2_Function_Vector_vtk(self):
     reference="hex_2D_o2_cell_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Function_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_Function_Vector.xml",reference)
  def test_hex_contact_2D_order2_Function_Tensor_vtk(self):
     reference="hex_2D_o2_cell_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Function_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_Function_Tensor.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o2_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o2_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o2_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o2_f_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o2_f_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o2_f_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("hex_contact_2D_order2_onFace_FunctionOnContactOne_Tensor.xml",reference)


  # ======================================================================================================================
  def test_hex_contact_3D_order1_ContinuousFunction_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_3D_order1_Solution_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Solution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_Solution_Scalar.xml",reference)
  def test_hex_contact_3D_order1_Solution_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Solution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_Solution_Vector.xml",reference)
  def test_hex_contact_3D_order1_Solution_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Solution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_Solution_Tensor.xml",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ReducedSolution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_3D_order1_Function_Scalar_vtk(self):
     reference="hex_3D_o1_cell_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Function_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_Function_Scalar.xml",reference)
  def test_hex_contact_3D_order1_Function_Vector_vtk(self):
     reference="hex_3D_o1_cell_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Function_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_Function_Vector.xml",reference)
  def test_hex_contact_3D_order1_Function_Tensor_vtk(self):
     reference="hex_3D_o1_cell_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Function_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_Function_Tensor.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o1_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o1_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o1_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o1_f_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o1_f_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o1_f_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order1_onFace_FunctionOnContactOne_Tensor.xml",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order2_ContinuousFunction_Scalar_vtk(self):
     reference="hex_3D_o2_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Vector_vtk(self):
     reference="hex_3D_o2_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Tensor_vtk(self):
     reference="hex_3D_o2_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_3D_order2_Solution_Scalar_vtk(self):
     reference="hex_3D_o2_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Solution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_Solution_Scalar.xml",reference)
  def test_hex_contact_3D_order2_Solution_Vector_vtk(self):
     reference="hex_3D_o2_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Solution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_Solution_Vector.xml",reference)
  def test_hex_contact_3D_order2_Solution_Tensor_vtk(self):
     reference="hex_3D_o2_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Solution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_Solution_Tensor.xml",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ReducedSolution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_3D_order2_Function_Scalar_vtk(self):
     reference="hex_3D_o2_cell_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Function_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_Function_Scalar.xml",reference)
  def test_hex_contact_3D_order2_Function_Vector_vtk(self):
     reference="hex_3D_o2_cell_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Function_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_Function_Vector.xml",reference)
  def test_hex_contact_3D_order2_Function_Tensor_vtk(self):
     reference="hex_3D_o2_cell_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Function_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_Function_Tensor.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o2_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o2_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o2_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o2_f_boundary_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o2_f_boundary_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o2_f_boundary_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("hex_contact_3D_order2_onFace_FunctionOnContactOne_Tensor.xml",reference)

class Test_DXFiles(Test_VisualizationInterface):
  # ======================================================================================================================
  def test_hex_2D_order2_dx(self):
     reference="hex_2D_o1.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     saveDX(FINLEY_WORKDIR_PATH+"hex_2D_order2.dx",domain=dom)
     self.check_dx("hex_2D_order2.dx",reference)

  def test_hex_2D_order2_AllPoints_Scalar_dx(self):
     reference="hex_2D_o1_node_3xs.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_2D_order2_AllPoints_Scalar.dx",data_r=x_r[0],data_n=x_n[0],data=x[0])
     self.check_dx("hex_2D_order2_AllPoints_Scalar.dx",reference)
  def test_hex_2D_order2_02Points_Scalar_dx(self):
     reference="hex_2D_o1_node_2xs.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Solution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_2D_order2_O2Points_Scalar.dx",data_n=x_n[0],data=x[0])
     self.check_dx("hex_2D_order2_O2Points_Scalar.dx",reference)
  def test_hex_2D_order2_2Cells_Scalar_dx(self):
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     try: 
        saveDX(FINLEY_WORKDIR_PATH+"hex_2D_order2_2Cells_Scalar.dx",data=x[0],data_b=x_b[0])
        self.fail("non-matching data not detected.")
     except StandardError:
        pass
  def test_hex_2D_order2_BoundrayPoint_Scalar_dx(self):
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     try: 
        saveDX(FINLEY_WORKDIR_PATH+"hex_2D_order2_BoundrayPoint_Scalar.dx",data=x[0],data_b=x_b[0])
        self.fail("non-matching data not detected.")
     except StandardError:
        pass
  def test_hex_2D_order2_Cells_AllData_dx(self):
     reference="hex_2D_o1_cell_all.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_2D_order2_Cells_AllData.dx",data_s=x[0],data_v=x[0]*[1.,2.],data_t=x[0]*[[11.,12.],[21.,22.]],data_t2=x[0]*[[-11.,-12.],[-21.,-22.]])
     self.check_dx("hex_2D_order2_Cells_AllData.dx",reference)

  def test_hex_2D_order2_CellsPoints_AllData_dx(self):
     reference="hex_2D_o1_cellnode_all.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_2D_order2.msh")
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_2D_order2_CellsPoints_AllData.dx",data_sp=x_p[0],
                                                     data_vp=x_p[0]*[1.,2.],
                                                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                                                     data_sc=x_c[0],
                                                     data_vc=x_c[0]*[1.,2.],
                                                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_2D_order2_CellsPoints_AllData.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order1_ContinuousFunction_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order1_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order1_Solution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Solution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order1_Solution_Scalar.dx",reference)
  def test_hex_contact_2D_order1_Solution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Solution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_Solution_Vector.dx",reference)
  def test_hex_contact_2D_order1_Solution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Solution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_Solution_Tensor.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order1_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ReducedSolution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_2D_order1_Function_Scalar_dx(self):
     reference="hex_2D_o1_cell_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Function_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order1_Function_Scalar.dx",reference)
  def test_hex_contact_2D_order1_Function_Vector_dx(self):
     reference="hex_2D_o1_cell_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Function_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_Function_Vector.dx",reference)
  def test_hex_contact_2D_order1_Function_Tensor_dx(self):
     reference="hex_2D_o1_cell_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_Function_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_Function_Tensor.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o1_boundary_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order1_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o1_boundary_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o1_boundary_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order1_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_FunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order2_ContinuousFunction_Scalar_dx(self):
     reference="hex_2D_o2_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order2_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Vector_dx(self):
     reference="hex_2D_o2_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Tensor_dx(self):
     reference="hex_2D_o2_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order2_Solution_Scalar_dx(self):
     reference="hex_2D_o2_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Solution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order2_Solution_Scalar.dx",reference)
  def test_hex_contact_2D_order2_Solution_Vector_dx(self):
     reference="hex_2D_o2_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Solution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_Solution_Vector.dx",reference)
  def test_hex_contact_2D_order2_Solution_Tensor_dx(self):
     reference="hex_2D_o2_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Solution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_Solution_Tensor.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order2_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ReducedSolution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_2D_order2_Function_Scalar_dx(self):
     reference="hex_2D_o2_cell_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Function_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order2_Function_Scalar.dx",reference)
  def test_hex_contact_2D_order2_Function_Vector_dx(self):
     reference="hex_2D_o2_cell_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Function_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_Function_Vector.dx",reference)
  def test_hex_contact_2D_order2_Function_Tensor_dx(self):
     reference="hex_2D_o2_cell_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_Function_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_Function_Tensor.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o2_boundary_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_2D_order2_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o2_boundary_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o2_boundary_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_2D_order2_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_FunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order1_ContinuousFunction_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order1_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order1_Solution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Solution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order1_Solution_Scalar.dx",reference)
  def test_hex_contact_3D_order1_Solution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Solution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_Solution_Vector.dx",reference)
  def test_hex_contact_3D_order1_Solution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Solution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_Solution_Tensor.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order1_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ReducedSolution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_3D_order1_Function_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Function_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order1_Function_Scalar.dx",reference)
  def test_hex_contact_3D_order1_Function_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Function_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_Function_Vector.dx",reference)
  def test_hex_contact_3D_order1_Function_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_Function_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_Function_Tensor.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order1_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order1_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_FunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order2_ContinuousFunction_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order2_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order2_Solution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Solution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order2_Solution_Scalar.dx",reference)
  def test_hex_contact_3D_order2_Solution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Solution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_Solution_Vector.dx",reference)
  def test_hex_contact_3D_order2_Solution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Solution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_Solution_Tensor.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order2_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ReducedSolution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_3D_order2_Function_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Function_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order2_Function_Scalar.dx",reference)
  def test_hex_contact_3D_order2_Function_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Function_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_Function_Vector.dx",reference)
  def test_hex_contact_3D_order2_Function_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_Function_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_Function_Tensor.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("hex_contact_3D_order2_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(FINLEY_TEST_MESH_PATH+"hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX(FINLEY_WORKDIR_PATH+"hex_contact_3D_order2_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_FunctionOnBoundary_Tensor.dx",reference)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_VTKFiles))
   suite.addTest(unittest.makeSuite(Test_DXFiles))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
