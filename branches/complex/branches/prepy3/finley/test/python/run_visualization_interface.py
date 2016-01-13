
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
from esys.escript import *
from esys.finley import ReadMesh

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")
# if os.name == "nt":
#    FINLEY_TEST_MESH_PATH = os.path.join(FINLEY_TEST_MESH_PATH,"win32")
FINLEY_WORKDIR_PATH=FINLEY_WORKDIR

class Test_VisualizationInterface(unittest.TestCase):
   def check_dx(self,f,reference_f):
      out_string=open(os.path.join(FINLEY_WORKDIR_PATH,f)).read().splitlines()
      ref_string=open(os.path.join(FINLEY_TEST_MESH_PATH,reference_f)).read().splitlines()
      c=0
      for l in range(0,len(ref_string)):
         if not ref_string[l].strip()[0]=="#":
           line=out_string[c].strip()
           if os.name == "nt":
               line=line.replace("e+00","e+0").replace("e-00","e-0")
           line=line.replace("e-00","e+00").replace("-0.000000e+00","0.000000e+00")
           self.assertEqual(line,ref_string[l].strip(),"line %d (%s) in dx file does not match reference (%s)"%(c,line,ref_string[l].strip()))
           c+=1

class Test_DXFiles(Test_VisualizationInterface):
  # ===========================================================================
  def test_hex_2D_order2_dx(self):
     reference="hex_2D_o1.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_2D_order2.msh"),optimize=False)
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_2D_order2.dx"),domain=dom)
     self.check_dx("hex_2D_order2.dx",reference)

  def test_hex_2D_order2_AllPoints_Scalar_dx(self):
     reference="hex_2D_o1_node_3xs.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_2D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_2D_order2_AllPoints_Scalar.dx"),data_r=x_r[0],data_n=x_n[0],data=x[0])
     self.check_dx("hex_2D_order2_AllPoints_Scalar.dx",reference)
  def test_hex_2D_order2_02Points_Scalar_dx(self):
     reference="hex_2D_o1_node_2xs.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_2D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_2D_order2_O2Points_Scalar.dx"),data_n=x_n[0],data=x[0])
     self.check_dx("hex_2D_order2_O2Points_Scalar.dx",reference)
  def test_hex_2D_order2_2Cells_Scalar_dx(self):
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     try: 
        saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_2D_order2_2Cells_Scalar.dx"),data=x[0],data_b=x_b[0])
        self.fail("non-matching data not detected.")
     except Exception:
        pass
  def test_hex_2D_order2_BoundaryPoint_Scalar_dx(self):
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_2D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     try: 
        saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_2D_order2_BoundaryPoint_Scalar.dx"),data=x[0],data_b=x_b[0])
        self.fail("non-matching data not detected.")
     except Exception:
        pass
  def test_hex_2D_order2_Cells_AllData_dx(self):
     reference="hex_2D_o1_cell_all.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_2D_order2_Cells_AllData.dx"),data_s=x[0],data_v=x[0]*[1.,2.],data_t=x[0]*[[11.,12.],[21.,22.]],data_t2=x[0]*[[-11.,-12.],[-21.,-22.]])
     self.check_dx("hex_2D_order2_Cells_AllData.dx",reference)

  def test_hex_2D_order2_CellsPoints_AllData_dx(self):
     reference="hex_2D_o1_cellnode_all.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_2D_order2.msh"),optimize=False)
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_2D_order2_CellsPoints_AllData.dx"),data_sp=x_p[0],
                                                     data_vp=x_p[0]*[1.,2.],
                                                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                                                     data_sc=x_c[0],
                                                     data_vc=x_c[0]*[1.,2.],
                                                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_2D_order2_CellsPoints_AllData.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order1_ContinuousFunction_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ContinuousFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order1_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ContinuousFunction_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_2D_order1_ContinuousFunction_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ContinuousFunction_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order1_Solution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_Solution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order1_Solution_Scalar.dx",reference)
  def test_hex_contact_2D_order1_Solution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_Solution_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_Solution_Vector.dx",reference)
  def test_hex_contact_2D_order1_Solution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_Solution_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_Solution_Tensor.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedSolution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order1_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedSolution_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_2D_order1_ReducedSolution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedSolution_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_2D_order1_Function_Scalar_dx(self):
     reference="hex_2D_o1_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_Function_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order1_Function_Scalar.dx",reference)
  def test_hex_contact_2D_order1_Function_Vector_dx(self):
     reference="hex_2D_o1_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_Function_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_Function_Vector.dx",reference)
  def test_hex_contact_2D_order1_Function_Tensor_dx(self):
     reference="hex_2D_o1_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_Function_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_Function_Tensor.dx",reference)
  def test_hex_contact_2D_order1_ReducedFunction_Scalar_dx(self):
     reference="hex_2D_o1_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order1_ReducedFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ReducedFunction_Vector_dx(self):
     reference="hex_2D_o1_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedFunction_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_ReducedFunction_Vector.dx",reference)
  def test_hex_contact_2D_order1_ReducedFunction_Tensor_dx(self):
     reference="hex_2D_o1_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedFunction_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_ReducedFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o1_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_FunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order1_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o1_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_FunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order1_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o1_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_FunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_FunctionOnBoundary_Tensor.dx",reference)
  def test_hex_contact_2D_order1_ReducedFunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o1_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedFunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order1_ReducedFunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order1_ReducedFunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o1_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedFunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order1_ReducedFunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order1_ReducedFunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o1_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order1_ReducedFunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order1_ReducedFunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_2D_order2_ContinuousFunction_Scalar_dx(self):
     reference="hex_2D_o2_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ContinuousFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order2_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Vector_dx(self):
     reference="hex_2D_o2_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ContinuousFunction_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_2D_order2_ContinuousFunction_Tensor_dx(self):
     reference="hex_2D_o2_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ContinuousFunction_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order2_Solution_Scalar_dx(self):
     reference="hex_2D_o2_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_Solution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order2_Solution_Scalar.dx",reference)
  def test_hex_contact_2D_order2_Solution_Vector_dx(self):
     reference="hex_2D_o2_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_Solution_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_Solution_Vector.dx",reference)
  def test_hex_contact_2D_order2_Solution_Tensor_dx(self):
     reference="hex_2D_o2_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_Solution_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_Solution_Tensor.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Scalar_dx(self):
     reference="hex_2D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedSolution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order2_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Vector_dx(self):
     reference="hex_2D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedSolution_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_2D_order2_ReducedSolution_Tensor_dx(self):
     reference="hex_2D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedSolution_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_2D_order2_Function_Scalar_dx(self):
     reference="hex_2D_o2_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_Function_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order2_Function_Scalar.dx",reference)
  def test_hex_contact_2D_order2_Function_Vector_dx(self):
     reference="hex_2D_o2_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_Function_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_Function_Vector.dx",reference)
  def test_hex_contact_2D_order2_Function_Tensor_dx(self):
     reference="hex_2D_o2_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_Function_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_Function_Tensor.dx",reference)
  def test_hex_contact_2D_order2_ReducedFunction_Scalar_dx(self):
     reference="hex_2D_o2_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order2_ReducedFunction_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ReducedFunction_Vector_dx(self):
     reference="hex_2D_o2_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedFunction_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_ReducedFunction_Vector.dx",reference)
  def test_hex_contact_2D_order2_ReducedFunction_Tensor_dx(self):
     reference="hex_2D_o2_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedFunction_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_ReducedFunction_Tensor.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o2_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_FunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order2_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o2_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_FunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order2_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o2_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_FunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_FunctionOnBoundary_Tensor.dx",reference)
  def test_hex_contact_2D_order2_ReducedFunctionOnBoundary_Scalar_dx(self):
     reference="hex_2D_o2_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedFunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_2D_order2_ReducedFunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_2D_order2_ReducedFunctionOnBoundary_Vector_dx(self):
     reference="hex_2D_o2_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedFunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.])
     self.check_dx("hex_contact_2D_order2_ReducedFunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_2D_order2_ReducedFunctionOnBoundary_Tensor_dx(self):
     reference="hex_2D_o2_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_2D_order2_ReducedFunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.],[21.,22.]])
     self.check_dx("hex_contact_2D_order2_ReducedFunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order1_ContinuousFunction_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ContinuousFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order1_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ContinuousFunction_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_3D_order1_ContinuousFunction_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ContinuousFunction_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order1_Solution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_Solution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order1_Solution_Scalar.dx",reference)
  def test_hex_contact_3D_order1_Solution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_Solution_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_Solution_Vector.dx",reference)
  def test_hex_contact_3D_order1_Solution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_Solution_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_Solution_Tensor.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedSolution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order1_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedSolution_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_3D_order1_ReducedSolution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedSolution_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_3D_order1_Function_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_Function_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order1_Function_Scalar.dx",reference)
  def test_hex_contact_3D_order1_Function_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_Function_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_Function_Vector.dx",reference)
  def test_hex_contact_3D_order1_Function_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_Function_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_Function_Tensor.dx",reference)
  def test_hex_contact_3D_order1_ReducedFunction_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order1_ReducedFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ReducedFunction_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedFunction_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_ReducedFunction_Vector.dx",reference)
  def test_hex_contact_3D_order1_ReducedFunction_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedFunction_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_ReducedFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_FunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order1_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_FunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order1_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_FunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_FunctionOnBoundary_Tensor.dx",reference)
  def test_hex_contact_3D_order1_ReducedFunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedFunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order1_ReducedFunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order1_ReducedFunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedFunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order1_ReducedFunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order1_ReducedFunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order1_ReducedFunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order1_ReducedFunctionOnBoundary_Tensor.dx",reference)
  # ======================================================================================================================
  def test_hex_contact_3D_order2_ContinuousFunction_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ContinuousFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order2_ContinuousFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ContinuousFunction_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_ContinuousFunction_Vector.dx",reference)
  def test_hex_contact_3D_order2_ContinuousFunction_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ContinuousFunction_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_ContinuousFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order2_Solution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_Solution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order2_Solution_Scalar.dx",reference)
  def test_hex_contact_3D_order2_Solution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_Solution_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_Solution_Vector.dx",reference)
  def test_hex_contact_3D_order2_Solution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_Solution_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_Solution_Tensor.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Scalar_dx(self):
     reference="hex_3D_o1_node_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedSolution_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order2_ReducedSolution_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Vector_dx(self):
     reference="hex_3D_o1_node_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedSolution_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_ReducedSolution_Vector.dx",reference)
  def test_hex_contact_3D_order2_ReducedSolution_Tensor_dx(self):
     reference="hex_3D_o1_node_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedSolution_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_ReducedSolution_Tensor.dx",reference)
  def test_hex_contact_3D_order2_Function_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_Function_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order2_Function_Scalar.dx",reference)
  def test_hex_contact_3D_order2_Function_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_Function_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_Function_Vector.dx",reference)
  def test_hex_contact_3D_order2_Function_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_Function_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_Function_Tensor.dx",reference)
  def test_hex_contact_3D_order2_ReducedFunction_Scalar_dx(self):
     reference="hex_3D_o1_cell_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedFunction_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order2_ReducedFunction_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ReducedFunction_Vector_dx(self):
     reference="hex_3D_o1_cell_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedFunction_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_ReducedFunction_Vector.dx",reference)
  def test_hex_contact_3D_order2_ReducedFunction_Tensor_dx(self):
     reference="hex_3D_o1_cell_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedFunction_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_ReducedFunction_Tensor.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_FunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order2_FunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_FunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_FunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order2_FunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_FunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_FunctionOnBoundary_Tensor.dx",reference)
  def test_hex_contact_3D_order2_ReducedFunctionOnBoundary_Scalar_dx(self):
     reference="hex_3D_o1_boundary_s.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedFunctionOnBoundary_Scalar.dx"),data=x[0])
     self.check_dx("hex_contact_3D_order2_ReducedFunctionOnBoundary_Scalar.dx",reference)
  def test_hex_contact_3D_order2_ReducedFunctionOnBoundary_Vector_dx(self):
     reference="hex_3D_o1_boundary_v.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedFunctionOnBoundary_Vector.dx"),data=x[0]*[1.,2.,3.])
     self.check_dx("hex_contact_3D_order2_ReducedFunctionOnBoundary_Vector.dx",reference)
  def test_hex_contact_3D_order2_ReducedFunctionOnBoundary_Tensor_dx(self):
     reference="hex_3D_o1_boundary_t.dx"
     dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     saveDX(os.path.join(FINLEY_WORKDIR_PATH,"hex_contact_3D_order2_ReducedFunctionOnBoundary_Tensor.dx"),data=x[0]*[[11.,12.,13.],[21.,22.,23],[31.,32.,33.]])
     self.check_dx("hex_contact_3D_order2_ReducedFunctionOnBoundary_Tensor.dx",reference)

if __name__ == '__main__':
   import sys
   suite = unittest.TestSuite()
   # saveDX is not MPI parallel
   if getMPISizeWorld() == 1: 
       suite.addTest(unittest.makeSuite(Test_DXFiles))
       pass
   else:
       print("Test_DXFiles is dropped as number of processors >1")
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)
