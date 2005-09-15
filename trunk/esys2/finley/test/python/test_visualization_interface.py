# $Id$

import sys
import unittest
from esys.escript import *
from esys.finley import ReadMesh

REERFENCE_FILE_DIR="test_meshes"

class Test_VisualizationInterface(unittest.TestCase):
   def check_vtk(self,f,reference_f):
      out_string=open(f).read().splitlines()
      ref_string=open(REERFENCE_FILE_DIR+"/"+reference_f).read().splitlines()
      c=0
      for l in range(0,len(ref_string)):
         if not ref_string[l].strip()[:2]=="<!":
           self.failUnlessEqual(out_string[c].strip().replace("e-00","e+00"),ref_string[l].strip(),"line %d (%s) in vtk files does not match reference (%s)"%(c,out_string[c].strip(),ref_string[l].strip()))
           c+=1

   def check_dx(self,f,reference_f):
      out_string=open(f).read().splitlines()
      ref_string=open(REERFENCE_FILE_DIR+"/"+reference_f).read().splitlines()
      c=0
      for l in range(0,len(ref_string)):
         if not ref_string[l].strip()[0]=="#":
           self.failUnlessEqual(out_string[c].strip().replace("e-00","e+00"),ref_string[l].strip(),"line %d (%s) in dx file does not match reference (%s)"%(c,out_string[c].strip(),ref_string[l].strip()))
           c+=1

class Test_VTKFiles(Test_VisualizationInterface):
  # ======================================================================================================================
  def test_hex_contact_2D_order1_ContinuousFunction_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_2D_order1_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_2D_order1_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_2D_order1_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_2D_order1_Solution_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_2D_order1_Solution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_Solution_Scalar.xml",reference)
  def test_hex_contact_2D_order1_Solution_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_2D_order1_Solution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_Solution_Vector.xml",reference)
  def test_hex_contact_2D_order1_Solution_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_2D_order1_Solution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_Solution_Tensor.xml",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_2D_order1_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_2D_order1_ReducedSolution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_2D_order1_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_2D_order1_Function_Scalar_vtk(self):
     reference="hex_2D_o1_cell_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_2D_order1_Function_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_Function_Scalar.xml",reference)
  def test_hex_contact_2D_order1_Function_Vector_vtk(self):
     reference="hex_2D_o1_cell_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_2D_order1_Function_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_Function_Vector.xml",reference)
  def test_hex_contact_2D_order1_Function_Tensor_vtk(self):
     reference="hex_2D_o1_cell_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_2D_order1_Function_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_Function_Tensor.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o1_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o1_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o1_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o1_f_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o1_f_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o1_f_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order1_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order1_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order1_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order1_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order1_onFace_FunctionOnContactOne_Tensor.xml",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order2_ContinuousFunction_Scalar_vtk(self):
     reference="hex_2D_o2_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_2D_order2_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Vector_vtk(self):
     reference="hex_2D_o2_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_2D_order2_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Tensor_vtk(self):
     reference="hex_2D_o2_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_2D_order2_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_2D_order2_Solution_Scalar_vtk(self):
     reference="hex_2D_o2_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_2D_order2_Solution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_Solution_Scalar.xml",reference)
  def test_hex_contact_2D_order2_Solution_Vector_vtk(self):
     reference="hex_2D_o2_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_2D_order2_Solution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_Solution_Vector.xml",reference)
  def test_hex_contact_2D_order2_Solution_Tensor_vtk(self):
     reference="hex_2D_o2_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_2D_order2_Solution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_Solution_Tensor.xml",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Scalar_vtk(self):
     reference="hex_2D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_2D_order2_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Vector_vtk(self):
     reference="hex_2D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_2D_order2_ReducedSolution_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Tensor_vtk(self):
     reference="hex_2D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_2D_order2_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_2D_order2_Function_Scalar_vtk(self):
     reference="hex_2D_o2_cell_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_2D_order2_Function_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_Function_Scalar.xml",reference)
  def test_hex_contact_2D_order2_Function_Vector_vtk(self):
     reference="hex_2D_o2_cell_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_2D_order2_Function_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_Function_Vector.xml",reference)
  def test_hex_contact_2D_order2_Function_Tensor_vtk(self):
     reference="hex_2D_o2_cell_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_2D_order2_Function_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_Function_Tensor.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o2_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o2_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o2_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_2D_o2_f_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_2D_o2_f_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_2D_o2_f_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order2_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order2_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_2D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_2D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_2D_order2_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_2D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_2D_order2_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_vtk("data.hex_contact_2D_order2_onFace_FunctionOnContactOne_Tensor.xml",reference)

  # ======================================================================================================================
  def test_hex_contact_3D_order1_ContinuousFunction_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_3D_order1_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_3D_order1_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_3D_order1_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_3D_order1_Solution_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_3D_order1_Solution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_Solution_Scalar.xml",reference)
  def test_hex_contact_3D_order1_Solution_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_3D_order1_Solution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_Solution_Vector.xml",reference)
  def test_hex_contact_3D_order1_Solution_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_3D_order1_Solution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_Solution_Tensor.xml",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_3D_order1_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_3D_order1_ReducedSolution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_3D_order1_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_3D_order1_Function_Scalar_vtk(self):
     reference="hex_3D_o1_cell_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_3D_order1_Function_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_Function_Scalar.xml",reference)
  def test_hex_contact_3D_order1_Function_Vector_vtk(self):
     reference="hex_3D_o1_cell_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_3D_order1_Function_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_Function_Vector.xml",reference)
  def test_hex_contact_3D_order1_Function_Tensor_vtk(self):
     reference="hex_3D_o1_cell_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_3D_order1_Function_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_Function_Tensor.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o1_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o1_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o1_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o1_f_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o1_f_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o1_f_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order1_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order1_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o1_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o1_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order1_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o1_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order1_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order1_onFace_FunctionOnContactOne_Tensor.xml",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order2_ContinuousFunction_Scalar_vtk(self):
     reference="hex_3D_o2_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_3D_order2_ContinuousFunction_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_ContinuousFunction_Scalar.xml",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Vector_vtk(self):
     reference="hex_3D_o2_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_3D_order2_ContinuousFunction_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_ContinuousFunction_Vector.xml",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Tensor_vtk(self):
     reference="hex_3D_o2_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveVTK("hex_contact_3D_order2_ContinuousFunction_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_ContinuousFunction_Tensor.xml",reference)
  def test_hex_contact_3D_order2_Solution_Scalar_vtk(self):
     reference="hex_3D_o2_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_3D_order2_Solution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_Solution_Scalar.xml",reference)
  def test_hex_contact_3D_order2_Solution_Vector_vtk(self):
     reference="hex_3D_o2_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_3D_order2_Solution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_Solution_Vector.xml",reference)
  def test_hex_contact_3D_order2_Solution_Tensor_vtk(self):
     reference="hex_3D_o2_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveVTK("hex_contact_3D_order2_Solution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_Solution_Tensor.xml",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Scalar_vtk(self):
     reference="hex_3D_o1_node_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_3D_order2_ReducedSolution_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_ReducedSolution_Scalar.xml",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Vector_vtk(self):
     reference="hex_3D_o1_node_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_3D_order2_ReducedSolution_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_ReducedSolution_Vector.xml",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Tensor_vtk(self):
     reference="hex_3D_o1_node_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveVTK("hex_contact_3D_order2_ReducedSolution_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_ReducedSolution_Tensor.xml",reference)
  def test_hex_contact_3D_order2_Function_Scalar_vtk(self):
     reference="hex_3D_o2_cell_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_3D_order2_Function_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_Function_Scalar.xml",reference)
  def test_hex_contact_3D_order2_Function_Vector_vtk(self):
     reference="hex_3D_o2_cell_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_3D_order2_Function_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_Function_Vector.xml",reference)
  def test_hex_contact_3D_order2_Function_Tensor_vtk(self):
     reference="hex_3D_o2_cell_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveVTK("hex_contact_3D_order2_Function_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_Function_Tensor.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o2_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o2_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o2_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnBoundary_Scalar_vtk(self):
     reference="hex_3D_o2_f_boundary_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnBoundary_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnBoundary_Scalar.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnBoundary_Vector_vtk(self):
     reference="hex_3D_o2_f_boundary_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnBoundary_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnBoundary_Vector.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnBoundary_Tensor_vtk(self):
     reference="hex_3D_o2_f_boundary_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnBoundary(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnBoundary_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnBoundary_Tensor.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactZero_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnContactZero_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnContactZero_Scalar.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactZero_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnContactZero_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnContactZero_Vector.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactZero_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactZero(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnContactZero_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnContactZero_Tensor.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order2_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order2_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_FunctionOnContactOne_Tensor.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactOne_Scalar_vtk(self):
     reference="hex_3D_o2_contact_s.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnContactOne_Scalar.xml",data=x[0])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnContactOne_Scalar.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactOne_Vector_vtk(self):
     reference="hex_3D_o2_contact_v.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnContactOne_Vector.xml",data=x[0]*[1.,2.,3.])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnContactOne_Vector.xml",reference)
  def test_hex_contact_3D_order2_onFace_FunctionOnContactOne_Tensor_vtk(self):
     reference="hex_3D_o2_contact_t.xml"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2_onFace.msh")
     x=FunctionOnContactOne(dom).getX()
     saveVTK("hex_contact_3D_order2_onFace_FunctionOnContactOne_Tensor.xml",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_vtk("data.hex_contact_3D_order2_onFace_FunctionOnContactOne_Tensor.xml",reference)

class Test_DXFiles(Test_VisualizationInterface):
  # ======================================================================================================================
  def test_hex_contact_2D_order1_ContinuousFunction_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_2D_order1_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order1_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_2D_order1_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order1_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_2D_order1_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order1_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order1_Solution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_2D_order1_Solution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order1_Solution_Scalar.dx",reference)
  def test_hex_contact_2D_order1_Solution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_2D_order1_Solution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order1_Solution_Vector.dx",reference)
  def test_hex_contact_2D_order1_Solution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_2D_order1_Solution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order1_Solution_Tensor.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_2D_order1_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order1_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_2D_order1_ReducedSolution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order1_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_2D_order1_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order1_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_2D_order1_Function_Scalar_dx(self):
     reference="hex_2D_o1_cell_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_2D_order1_Function_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order1_Function_Scalar.dx",reference)
  def test_hex_contact_2D_order1_Function_Vector_dx(self):
     reference="hex_2D_o1_cell_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_2D_order1_Function_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order1_Function_Vector.dx",reference)
  def test_hex_contact_2D_order1_Function_Tensor_dx(self):
     reference="hex_2D_o1_cell_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_2D_order1_Function_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order1_Function_Tensor.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o1_boundary_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_2D_order1_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order1_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o1_boundary_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_2D_order1_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order1_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o1_boundary_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_2D_order1_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order1_FunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order2_ContinuousFunction_Scalar_dx(self):
     reference="hex_2D_o2_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_2D_order2_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order2_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Vector_dx(self):
     reference="hex_2D_o2_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_2D_order2_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order2_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Tensor_dx(self):
     reference="hex_2D_o2_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_2D_order2_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order2_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order2_Solution_Scalar_dx(self):
     reference="hex_2D_o2_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_2D_order2_Solution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order2_Solution_Scalar.dx",reference)
  def test_hex_contact_2D_order2_Solution_Vector_dx(self):
     reference="hex_2D_o2_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_2D_order2_Solution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order2_Solution_Vector.dx",reference)
  def test_hex_contact_2D_order2_Solution_Tensor_dx(self):
     reference="hex_2D_o2_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_2D_order2_Solution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order2_Solution_Tensor.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_2D_order2_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order2_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_2D_order2_ReducedSolution_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order2_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_2D_order2_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order2_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_2D_order2_Function_Scalar_dx(self):
     reference="hex_2D_o2_cell_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_2D_order2_Function_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order2_Function_Scalar.dx",reference)
  def test_hex_contact_2D_order2_Function_Vector_dx(self):
     reference="hex_2D_o2_cell_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_2D_order2_Function_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order2_Function_Vector.dx",reference)
  def test_hex_contact_2D_order2_Function_Tensor_dx(self):
     reference="hex_2D_o2_cell_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_2D_order2_Function_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order2_Function_Tensor.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o2_boundary_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_2D_order2_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_2D_order2_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o2_boundary_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_2D_order2_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.])
     self.check_dx("data.hex_contact_2D_order2_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o2_boundary_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_2D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_2D_order2_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("data.hex_contact_2D_order2_FunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order1_ContinuousFunction_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_3D_order1_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order1_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_3D_order1_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order1_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_3D_order1_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order1_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order1_Solution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_3D_order1_Solution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order1_Solution_Scalar.dx",reference)
  def test_hex_contact_3D_order1_Solution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_3D_order1_Solution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order1_Solution_Vector.dx",reference)
  def test_hex_contact_3D_order1_Solution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_3D_order1_Solution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order1_Solution_Tensor.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_3D_order1_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order1_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_3D_order1_ReducedSolution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order1_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_3D_order1_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order1_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_3D_order1_Function_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_3D_order1_Function_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order1_Function_Scalar.dx",reference)
  def test_hex_contact_3D_order1_Function_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_3D_order1_Function_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order1_Function_Vector.dx",reference)
  def test_hex_contact_3D_order1_Function_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_3D_order1_Function_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order1_Function_Tensor.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_3D_order1_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order1_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_3D_order1_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order1_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order1.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_3D_order1_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order1_FunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order2_ContinuousFunction_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_3D_order2_ContinuousFunction_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order2_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_3D_order2_ContinuousFunction_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order2_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ContinuousFunction(dom).getX()
     saveDX("hex_contact_3D_order2_ContinuousFunction_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order2_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order2_Solution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_3D_order2_Solution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order2_Solution_Scalar.dx",reference)
  def test_hex_contact_3D_order2_Solution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_3D_order2_Solution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order2_Solution_Vector.dx",reference)
  def test_hex_contact_3D_order2_Solution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Solution(dom).getX()
     saveDX("hex_contact_3D_order2_Solution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order2_Solution_Tensor.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_3D_order2_ReducedSolution_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order2_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_3D_order2_ReducedSolution_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order2_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=ReducedSolution(dom).getX()
     saveDX("hex_contact_3D_order2_ReducedSolution_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order2_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_3D_order2_Function_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_3D_order2_Function_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order2_Function_Scalar.dx",reference)
  def test_hex_contact_3D_order2_Function_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_3D_order2_Function_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order2_Function_Vector.dx",reference)
  def test_hex_contact_3D_order2_Function_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=Function(dom).getX()
     saveDX("hex_contact_3D_order2_Function_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order2_Function_Tensor.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_3D_order2_FunctionOnBoundary_Scalar.dx",data=x[0])
     self.check_dx("data.hex_contact_3D_order2_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_3D_order2_FunctionOnBoundary_Vector.dx",data=x[0]*[1.,2.,3.])
     self.check_dx("data.hex_contact_3D_order2_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(REERFENCE_FILE_DIR+"/hex_contact_3D_order2.msh")
     x=FunctionOnBoundary(dom).getX()
     saveDX("hex_contact_3D_order2_FunctionOnBoundary_Tensor.dx",data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("data.hex_contact_3D_order2_FunctionOnBoundary_Tensor.dx",reference)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_VTKFiles))
   suite.addTest(unittest.makeSuite(Test_DXFiles))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
