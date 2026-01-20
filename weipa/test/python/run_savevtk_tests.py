
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os, math
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from xml.dom import minidom
from esys.escript import ContinuousFunction, Function, ReducedFunction,\
            FunctionOnBoundary, ReducedFunctionOnBoundary,\
            FunctionOnContactZero, ReducedFunctionOnContactZero,\
            FunctionOnContactOne, ReducedFunctionOnContactOne,\
            Solution, ReducedSolution, getMPISizeWorld
# from esys.weipa import saveVTK

try:
    from esys import finley
    finleyInstalled = True
except:
    finleyInstalled = False
try:
    from esys import ripley
    ripleyInstalled = True
except:
    ripleyInstalled = False
try:
    from esys import oxley
    oxleyInstalled = True
except:
    oxleyInstalled = False

try:
     WEIPA_TEST_MESHES=os.environ['WEIPA_TEST_MESHES']
except KeyError:
     WEIPA_TEST_MESHES='meshes'

try:
     WEIPA_WORKDIR=os.environ['WEIPA_WORKDIR']
except KeyError:
     WEIPA_WORKDIR='.'

class VTKParser():
    """
    VTU (VTK XML format) file parser.
    """
    def __init__(self):
        self.doc=None

    def parse(self, filename):
        # load file and remove superfluous whitespace
        dom=minidom.parse(filename)
        dom=minidom.parseString(dom.toxml()
                .replace('>\n','>').replace('\n<','<').replace('\n',' '))
        self.doc = dom.documentElement
        # check for general file validity
        return self.doc.tagName=='VTKFile' and \
               self.doc.getAttribute('type')=='UnstructuredGrid'

    def getTimeAndCycle(self):
        time=None
        cycle=None
        for e in self.doc.getElementsByTagName('FieldData'):
            for d in e.childNodes:
                if d.getAttribute('Name')=='TIME':
                    time = float(d.childNodes[0].data)
                elif d.getAttribute('Name')=='CYCLE':
                    cycle = int(d.childNodes[0].data)
        return time, cycle

    def getNumPointsAndCells(self):
        piece = self.doc.getElementsByTagName('Piece')[0]
        return int(piece.getAttribute('NumberOfPoints')), int(piece.getAttribute('NumberOfCells'))

    def getPoints(self):
        points = self.doc.getElementsByTagName('Points')[0]
        array = points.childNodes[0]
        nComp = int(array.getAttribute('NumberOfComponents'))
        flatlist = list(map(float, array.childNodes[0].data.split()))
        return [flatlist[i:i+nComp] for i in range(0,len(flatlist),nComp)]

    def getCellConnectivity(self):
        conn=None
        cells = self.doc.getElementsByTagName('Cells')[0]
        for d in cells.childNodes:
            if d.getAttribute('Name')=='connectivity':
                conn = list(map(int, d.childNodes[0].data.split()))
        return conn

    def getCellOffsets(self):
        offsets=None
        cells = self.doc.getElementsByTagName('Cells')[0]
        for d in cells.childNodes:
            if d.getAttribute('Name')=='offsets':
                offsets = list(map(int, d.childNodes[0].data.split()))
        return offsets

    def getCellTypes(self):
        types=None
        cells = self.doc.getElementsByTagName('Cells')[0]
        for d in cells.childNodes:
            if d.getAttribute('Name')=='types':
                types = list(map(int, d.childNodes[0].data.split()))
        return types

    def getPointData(self):
        data={}
        pdata = self.doc.getElementsByTagName('PointData')
        if len(pdata)>0:
            for d in pdata[0].childNodes:
                if d.hasChildNodes():
                    nComp = int(d.getAttribute('NumberOfComponents'))
                    flatlist = list(map(float, d.childNodes[0].data.split()))
                    if nComp==1:
                        data[d.getAttribute('Name')] = flatlist
                    else:
                        l=[flatlist[i:i+nComp] for i in range(0,len(flatlist),nComp)]
                        data[d.getAttribute('Name')] = l
        return data

    def getCellData(self):
        data={}
        cdata = self.doc.getElementsByTagName('CellData')
        if len(cdata)>0:
            for d in cdata[0].childNodes:
                if d.hasChildNodes():
                    nComp = int(d.getAttribute('NumberOfComponents'))
                    flatlist = list(map(float, d.childNodes[0].data.split()))
                    if nComp==1:
                        data[d.getAttribute('Name')] = flatlist
                    else:
                        l=[flatlist[i:i+nComp] for i in range(0,len(flatlist),nComp)]
                        data[d.getAttribute('Name')] = l
        return data


# class Test_VTKSaver(unittest.TestCase):
#     # Compares two lists of numbers using the L2 norm, returns true if they
#     # match up to a tolerance TOL, false otherwise
#     def numericCompareL2(self, vector1, vector2):
#         if vector2 == None: return False
#         TOL = 2.0e-5
#         if len(vector1) != len(vector2): return False
#         diff = 0.0
#         for i in range(0, len(vector1)):
#             tmp = vector1[i] - vector2[i]
#             diff += tmp * tmp
#         if math.sqrt(diff) > TOL: return False
#         return True

#     # Compares two data arrays (of scalars, vectors & tensors) which may be
#     # permuted. IndexMap is used to map one value in d1 to one or more values
#     # in d2.
#     def compareDataWithMap(self, d1, d2, indexMap):
#         if len(d1) != len(d2): return False
#         for i in range(len(d1)):
#             if i in indexMap:
#                 jlist=indexMap[i]
#                 if type(d1[i])==list:
#                     if not max([ self.numericCompareL2(d1[i], d2[j]) for j in jlist ]):
#                         return False
#                 else:
#                     if not max([ self.numericCompareL2([d1[i]], [d2[j]]) for j in jlist ]):
#                         return False
#             else:
#                 return False
#         return True

#     # Compares two VTK files, asserts if any element is not equal.
#     # file2 is the reference file to compare against
#     def compareVTKfiles(self, file1, file2):
#         p1=VTKParser()
#         p2=VTKParser()
#         self.assertTrue(p1.parse(file1), "Invalid vtu file")
#         p2.parse(file2)
#         self.assertEqual(p1.getTimeAndCycle(), p2.getTimeAndCycle())
#         nPoints1, nCells1 = p1.getNumPointsAndCells()
#         self.assertEqual((nPoints1, nCells1), p2.getNumPointsAndCells())

#         nodeList1 = p1.getPoints()
#         nodeList2 = p2.getPoints()
#         self.assertEqual(len(nodeList1), nPoints1)

#         # Find mapping of nodes in file 1 to file 2 (they may be permuted)
#         nodeMap1to2 = {}
#         for i1 in range(len(nodeList1)):
#             indexList=[]
#             for i2 in range(len(nodeList2)):
#                 if self.numericCompareL2(nodeList1[i1], nodeList2[i2]):
#                     indexList.append(i2)
#             self.assertNotEqual(len(indexList), 0,
#                              "Node with coordinates %s missing."%nodeList1[i1])
#             nodeMap1to2[i1]=indexList

#         # cells
#         offsets1 = p1.getCellOffsets()
#         offsets2 = p2.getCellOffsets()
#         self.assertEqual(len(offsets1), nCells1)
#         types1 = p1.getCellTypes()
#         types2 = p2.getCellTypes()
#         self.assertEqual(len(types1), nCells1)
#         conn1 = p1.getCellConnectivity()
#         conn2 = p2.getCellConnectivity()
#         elementList1=[]
#         elementList2=[]
#         lastOffset1=0
#         lastOffset2=0
#         for i in range(nCells1):
#             elementList1.append(conn1[lastOffset1:offsets1[i]])
#             elementList2.append(conn2[lastOffset2:offsets2[i]])
#             lastOffset1=offsets1[i]
#             lastOffset2=offsets2[i]
#         self.assertEqual(len(elementList1), len(elementList2))

#         # find element mapping, then compare types.
#         # offsets are compared implicitly by creating the mapping
#         elementMap1to2 = {}
#         for i1 in range(len(elementList1)):
#             index=None
#             for i2 in range(len(elementList2)):
#                 mappedL1=[nodeMap1to2[i] for i in elementList1[i1]]
#                 if len(mappedL1)==len(elementList2[i2]):
#                     matches=[elementList2[i2][i] in mappedL1[i] for i in range(len(elementList2[i2]))]
#                     if min(matches):
#                         index=i2
#                         break
#             self.assertNotEqual(index, None,
#                                  "Element %s is missing."%elementList1[i1])
#             elementMap1to2[i1]=[index]

#         types1n = [types1[elementMap1to2[i][0]] for i in range(nCells1)]
#         self.assertEqual(types1n, types2)

#         # point data
#         pdata1=p1.getPointData()
#         pdata2=p2.getPointData()
#         self.assertEqual(len(pdata1), len(pdata2))
#         for name in pdata2:
#             self.assertTrue(name in pdata1, "Point variable '%s' missing"%name)
#             self.assertEqual(len(pdata1[name]), nPoints1)
#             if not name.startswith('mesh_vars/'):
#                 self.assertTrue(self.compareDataWithMap(
#                     pdata1[name], pdata2[name], nodeMap1to2),
#                     "Point data in '%s' does not match" % name)

#         # cell data
#         cdata1=p1.getCellData()
#         cdata2=p2.getCellData()
#         self.assertEqual(len(cdata1), len(cdata2))
#         for name in cdata2:
#             self.assertTrue(name in cdata1, "Cell variable '%s' missing"%name)
#             self.assertEqual(len(cdata1[name]), nCells1)
#             if not name.startswith('mesh_vars/'):
#                 self.assertTrue(self.compareDataWithMap(
#                     cdata1[name], cdata2[name], elementMap1to2),
#                     "Cell data in '%s' does not match" % name)

#     def check_vtk(self, reference, fspaces=[], **data):
#         outFileBase="out_"+reference
#         saveVTK(os.path.join(WEIPA_WORKDIR, outFileBase), write_meshdata=True, **data)
#         if len(fspaces)>0:
#             for fs in fspaces:
#                 ref=os.path.join(WEIPA_TEST_MESHES, reference+"_"+fs+".vtu")
#                 out=os.path.join(WEIPA_WORKDIR, outFileBase+"_"+fs+".vtu")
#                 self.compareVTKfiles(out, ref)
#         else:
#             ref=os.path.join(WEIPA_TEST_MESHES, reference+".vtu")
#             out=os.path.join(WEIPA_WORKDIR, outFileBase+".vtu")
#             self.compareVTKfiles(out, ref)


# @unittest.skipIf(not finleyInstalled, "Skipping finley saveVTK tests since finley not installed")
# class Test_Finley_SaveVTK(Test_VTKSaver):

#   # === METADATA =============================================================

#   def test_metadata_full(self):
#      fn=os.path.join(WEIPA_WORKDIR, "metadata0.vtu")
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES, "hex_2D_order2.msh"), optimize=False)
#      saveVTK(fn, x=dom.getX(), metadata_schema={"gml":"http://www.opengis.net/gml"}, metadata='<dummy>hello world</dummy><timeStamp uom="s">1234</timeStamp>')
#      # testing:
#      dom=VTKParser()
#      self.assertTrue(dom.parse(fn), "Invalid vtu file")
#      self.assertEqual(dom.doc.getAttribute("xmlns:gml"), "http://www.opengis.net/gml")
#      self.assertTrue(len(dom.doc.getElementsByTagName("MetaData"))>0, "No MetaData element found")
#      mdNode = dom.doc.getElementsByTagName("MetaData")[0]
#      dummy = None
#      timestamp = None
#      for n in mdNode.childNodes:
#          if n.tagName == 'dummy':
#              dummy = n.childNodes[0].data
#              self.assertEqual(dummy, 'hello world')
#          elif n.tagName == 'timeStamp':
#              uom = n.getAttribute('uom')
#              self.assertEqual(uom, 's')
#              timestamp = n.childNodes[0].data
#              self.assertEqual(timestamp, '1234')
#          else:
#              self.fail('Unexpected metadata tag found')
#      self.assertNotEqual(dummy, None)
#      self.assertNotEqual(timestamp, None)

#   def test_metadata_no_schema(self):
#      fn=os.path.join(WEIPA_WORKDIR, "metadata1.vtu")
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES, "hex_2D_order2.msh"), optimize=False)
#      saveVTK(fn, x=dom.getX(), metadata='<dummy>hello world</dummy><timeStamp uom="s">1234</timeStamp>')
#      # testing:
#      dom=VTKParser()
#      self.assertTrue(dom.parse(fn), "Invalid vtu file")
#      self.assertTrue(len(dom.doc.getElementsByTagName("MetaData"))>0, "No MetaData element found")
#      mdNode = dom.doc.getElementsByTagName("MetaData")[0]
#      dummy=None
#      timeStamp=None
#      for n in mdNode.childNodes:
#          if n.tagName == 'dummy':
#              dummy = n.childNodes[0].data
#              self.assertEqual(dummy, 'hello world')
#          elif n.tagName == 'timeStamp':
#              uom = n.getAttribute('uom')
#              self.assertEqual(uom, 's')
#              timestamp = n.childNodes[0].data
#              self.assertEqual(timestamp, '1234')
#          else:
#              self.fail('Unexpected metadata tag found')
#      self.assertNotEqual(dummy, None)
#      self.assertNotEqual(timestamp, None)

#   def test_metadata_schema_only(self):
#      fn=os.path.join(WEIPA_WORKDIR, "metadata2.vtu")
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES, "hex_2D_order2.msh"), optimize=False)
#      saveVTK(fn, x=dom.getX(), metadata_schema={"gml":"http://www.opengis.net/gml"})
#      # testing:
#      dom=VTKParser()
#      self.assertTrue(dom.parse(fn), "Invalid vtu file")
#      self.assertEqual(dom.doc.getAttribute("xmlns:gml"), "http://www.opengis.net/gml")

#   # === Finley hex 2D order 1 with contacts ===================================

#   def test_hex_contact_2D_order1_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_onFace_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                             data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_onFace_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                             data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_onFace_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_onFace_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_onFace_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order1_onFace_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Finley hex 2D order 2 with contacts ===================================

#   def test_hex_contact_2D_order2_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_2D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_2D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_2D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_onFace_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                             data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_onFace_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                             data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o2_rcontact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_onFace_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                            data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_onFace_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                            data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o2_rcontact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_onFace_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                            data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_contact_2D_order2_onFace_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
#                                            data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Finley hex 2D order 2 =================================================

#   def test_hex_2D_order2_empty(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
#      self.check_vtk("hex_2D_o2", domain=dom)

#   def test_hex_2D_order2_AllPoints_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
#      x=Solution(dom).getX()
#      x_r=ReducedSolution(dom).getX()
#      x_n=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o2_node_3xs", ['Elements','ReducedElements'], data_r=x_r[0], data_n=x_n[0], data=x[0])

#   def test_hex_2D_order2_2Cells_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
#      x=Function(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2_cell_2xs", ['ReducedElements','ReducedFaceElements'], data_b=x_b[0], data=x[0])

#   def test_hex_2D_order2_BoundaryPoint_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2_boundary_2xs", ['Elements','ReducedFaceElements'], data=x[0],data_b=x_b[0])

#   def test_hex_2D_order2_Cells_AllData(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_2D_o2_cell_all",
#                     data_s=x[0],
#                     data_v=x[0]*[1.,2.],
#                     data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2_CellsPoints_AllData(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
#      x_c=Function(dom).getX()
#      x_p=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o2_cellnode_all", ['Elements','ReducedElements'],
#                     data_sp=x_p[0],
#                     data_vp=x_p[0]*[1.,2.],
#                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
#                     data_sc=x_c[0],
#                     data_vc=x_c[0]*[1.,2.],
#                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

#   # === Finley hex 2D order 2 (full) ==========================================

#   def test_hex_2D_order2p_empty(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      self.check_vtk("hex_2D_o2p", domain=dom)

#   def test_hex_2D_order2p_AllPoints_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=Solution(dom).getX()
#      x_r=ReducedSolution(dom).getX()
#      x_n=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o1_node_3xs", ['Elements','ReducedElements'], data_r=x_r[0], data_n=x_n[0], data=x[0])

#   def test_hex_2D_order2p_2Cells_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=Function(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2p_cell_2xs", ['Elements','ReducedFaceElements'], data_b=x_b[0], data=x[0])

#   def test_hex_2D_order2p_BoundaryPoint_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2p_boundary_2xs", ['Elements','ReducedFaceElements'], data=x[0],data_b=x_b[0])

#   def test_hex_2D_order2p_CellsPoints_AllData(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x_c=Function(dom).getX()
#      x_p=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o2p_cellnode_all",
#                     data_sp=x_p[0],
#                     data_vp=x_p[0]*[1.,2.],
#                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
#                     data_sc=x_c[0],
#                     data_vc=x_c[0]*[1.,2.],
#                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2p_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2p_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2p_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_2D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
#                                         data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2p_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_2D_o2p_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2p_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_2D_o2p_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                         data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2p_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2p_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                            data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_order2p_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_o2p_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                            data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Finley hex 2D macro ===================================================

#   def test_hex_2D_macro_empty(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      self.check_vtk("hex_2D_o2p", domain=dom)

#   def test_hex_2D_macro_AllPoints_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=Solution(dom).getX()
#      x_r=ReducedSolution(dom).getX()
#      x_n=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o1_node_3xs", ['Elements','ReducedElements'], data_r=x_r[0], data_n=x_n[0], data=x[0])

#   def test_hex_2D_macro_CellsPoints(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x_c=Function(dom).getX()
#      x_p=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_macro_cellnode_all",
#                     data_sp=x_p[0],
#                     data_vp=x_p[0]*[1.,2.],
#                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
#                     data_sc=x_c[0],
#                     data_vc=x_c[0]*[1.,2.],
#                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_macro_2Cells_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=Function(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_cell_2xs", ['Elements','FaceElements'], data_b=x_b[0], data=x[0])

#   def test_hex_2D_macro_BoundaryPoint_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_boundary_2xs", ['Elements','FaceElements'], data_b=x_b[0], data=x[0])

#   def test_hex_2D_macro_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_macro_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_macro_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_2D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
#                                         data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_macro_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_2D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_macro_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_2D_macro_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_macro_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                              data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_hex_2D_macro_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                              data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Finley hex 3D order 1 with contacts ===================================

#   def test_hex_contact_3D_order1_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_onFace_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_onFace_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_onFace_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_onFace_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_onFace_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order1_onFace_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   # === Finley hex 3D order 2 with contacts ===================================

#   def test_hex_contact_3D_order2_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_3D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_3D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_3D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_onFace_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_onFace_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_onFace_FunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
#      x=FunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_onFace_ReducedFunctionOnContactZero(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactZero(dom).getX()
#      self.check_vtk("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_onFace_FunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
#      x=FunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_contact_3D_order2_onFace_ReducedFunctionOnContactOne(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
#      x=ReducedFunctionOnContactOne(dom).getX()
#      self.check_vtk("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   # === Finley hex 3D order 2 (full) ==========================================

#   def test_hex_3D_order2p_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_order2p_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_order2p_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_3D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                         data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_order2p_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_3D_o2p_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_order2p_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_3D_o2p_rcell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                         data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_order2p_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o2p_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_order2p_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_o2p_rboundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   # === Finley hex 3D macro ===================================================

#   def test_hex_3D_macro_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_macro_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_macro_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("hex_3D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                         data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_macro_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("hex_3D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_macro_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("hex_3D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_macro_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                              data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_hex_3D_macro_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("hex_3D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                              data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   # === Finley tet 2D order 1 =================================================

#   def test_tet_2D_order1_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order1_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("tet_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order1_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("tet_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order1_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("tet_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order1_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("tet_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order1_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order1_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Finley tet 2D order 2 =================================================

#   def test_tet_2D_order2(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      self.check_vtk("tet_2D_o2", domain=dom)

#   def test_tet_2D_order2_AllPoints_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=Solution(dom).getX()
#      x_r=ReducedSolution(dom).getX()
#      x_n=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o1_node_3xs", ['Elements','ReducedElements'], data_r=x_r[0], data_n=x_n[0], data=x[0])

#   def test_tet_2D_order2_02Points_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=Solution(dom).getX()
#      x_n=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o2_node_2xs", data_n=x_n[0], data=x[0])

#   def test_tet_2D_order2_2Cells_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=Function(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_o2_cell_2xs", ['Elements','ReducedFaceElements'], data_b=x_b[0], data=x[0])

#   def test_tet_2D_order2_BoundaryPoint_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_o2_boundary_2xs", ['Elements','ReducedFaceElements'], data=x[0],data_b=x_b[0])

#   def test_tet_2D_order2_Cells_AllData(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("tet_2D_o2_cell_all",
#                     data_s=x[0],
#                     data_v=x[0]*[1.,2.],
#                     data_t=x[0]*[[11.,12.],[21.,22.]],
#                     data_t2=x[0]*[[-11.,-12.],[-21.,-22.]])

#   def test_tet_2D_order2_CellsPoints_AllData(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x_c=Function(dom).getX()
#      x_p=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o2_cellnode_all",
#                     data_sp=x_p[0],
#                     data_vp=x_p[0]*[1.,2.],
#                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
#                     data_sc=x_c[0],
#                     data_vc=x_c[0]*[1.,2.],
#                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order2_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order2_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order2_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("tet_2D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order2_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("tet_2D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order2_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("tet_2D_o2_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order2_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_order2_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Finley tet 2D macro ===================================================

#   def test_tet_2D_macro(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      self.check_vtk("tet_2D_o2", domain=dom)

#   def test_tet_2D_macro_AllPoints_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=Solution(dom).getX()
#      x_r=ReducedSolution(dom).getX()
#      x_n=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o1_node_3xs", ['Elements','ReducedElements'], data_r=x_r[0], data_n=x_n[0], data=x[0])

#   def test_tet_2D_macro_02Points_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=Solution(dom).getX()
#      x_n=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o2_node_2xs", data_n=x_n[0], data=x[0])

#   def test_tet_2D_macro_2Cells_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=Function(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_cell_2xs", ['Elements','FaceElements'], data_b=x_b[0], data=x[0])

#   def test_tet_2D_macro_BoundaryPoint_Scalar(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      x_b=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_boundary_2xs", ['Elements','FaceElements'], data_b=x_b[0], data=x[0])

#   def test_tet_2D_macro_CellsPoints_AllData(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x_c=Function(dom).getX()
#      x_p=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_macro_cellnode_all",
#                     data_sp=x_p[0],
#                     data_vp=x_p[0]*[1.,2.],
#                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
#                     data_sc=x_c[0],
#                     data_vc=x_c[0]*[1.,2.],
#                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_macro_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_macro_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_macro_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("tet_2D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
#                                        data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_macro_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("tet_2D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                          data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_macro_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("tet_2D_macro_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_macro_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                              data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_tet_2D_macro_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                              data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Finley tet 3D order 1 =================================================

#   def test_tet_3D_order1_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order1_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("tet_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order1_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("tet_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order1_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("tet_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order1_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("tet_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order1_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order1_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   # === Finley tet 3D order 2 =================================================

#   def test_tet_3D_order2_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order2_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order2_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("tet_3D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order2_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("tet_3D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order2_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("tet_3D_o2_rcell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order2_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_3D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_order2_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_3D_o2_rboundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   # === Finley tet 3D macro ===================================================

#   def test_tet_3D_macro_ContinuousFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_macro_Solution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
#      x=Solution(dom).getX()
#      self.check_vtk("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_macro_ReducedSolution(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("tet_3D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_macro_Function(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
#      x=Function(dom).getX()
#      self.check_vtk("tet_3D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_macro_ReducedFunction(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("tet_3D_macro_rcell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_macro_FunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_3D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                              data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_tet_3D_macro_ReducedFunctionOnBoundary(self):
#      dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("tet_3D_macro_rboundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                               data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])



# @unittest.skipIf(getMPISizeWorld()>4, "Skipping ripley saveVTK tests since MPI size > 4")
# @unittest.skipIf(not ripleyInstalled, "Skipping ripley saveVTK tests since ripley not installed")
# class Test_Ripley_SaveVTK(Test_VTKSaver):

#   # === Ripley 2D =============================================================

#   def test_ripley_2D_ContinuousFunction(self):
#      dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("ripley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_ripley_2D_Solution(self):
#      dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=Solution(dom).getX()
#      self.check_vtk("ripley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_ripley_2D_ReducedSolution(self):
#      dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("ripley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_ripley_2D_Function(self):
#      dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=Function(dom).getX()
#      self.check_vtk("ripley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_ripley_2D_ReducedFunction(self):
#      dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("ripley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_ripley_2D_FunctionOnBoundary(self):
#      dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("ripley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_ripley_2D_ReducedFunctionOnBoundary(self):
#      dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("ripley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === Ripley 3D =============================================================

#   def test_ripley_3D_ContinuousFunction(self):
#      dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("ripley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_ripley_3D_Solution(self):
#      dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=Solution(dom).getX()
#      self.check_vtk("ripley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_ripley_3D_ReducedSolution(self):
#      dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("ripley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_ripley_3D_Function(self):
#      dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=Function(dom).getX()
#      self.check_vtk("ripley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_ripley_3D_ReducedFunction(self):
#      dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("ripley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_ripley_3D_FunctionOnBoundary(self):
#      dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("ripley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_ripley_3D_ReducedFunctionOnBoundary(self):
#      dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("ripley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

# @unittest.skipIf(getMPISizeWorld()>4, "Skipping oxley saveVTK tests since MPI size > 4")
# @unittest.skipIf(not oxleyInstalled, "Skipping oxley saveVTK tests since oxley not installed")
# class Test_Oxley_SaveVTK(Test_VTKSaver):

#   # === oxley 2D =============================================================

#   def test_oxley_2D_ContinuousFunction(self):
#      dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("oxley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_oxley_2D_Solution(self):
#      dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=Solution(dom).getX()
#      self.check_vtk("oxley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_oxley_2D_ReducedSolution(self):
#      dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("oxley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_oxley_2D_Function(self):
#      dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=Function(dom).getX()
#      self.check_vtk("oxley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_oxley_2D_ReducedFunction(self):
#      dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("oxley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
#                                       data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_oxley_2D_FunctionOnBoundary(self):
#      dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("oxley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   def test_oxley_2D_ReducedFunctionOnBoundary(self):
#      dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("oxley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
#                                           data_t=x[0]*[[11.,12.],[21.,22.]])

#   # === oxley 3D =============================================================

#   def test_oxley_3D_ContinuousFunction(self):
#      dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ContinuousFunction(dom).getX()
#      self.check_vtk("oxley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_oxley_3D_Solution(self):
#      dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=Solution(dom).getX()
#      self.check_vtk("oxley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_oxley_3D_ReducedSolution(self):
#      dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ReducedSolution(dom).getX()
#      self.check_vtk("oxley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_oxley_3D_Function(self):
#      dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=Function(dom).getX()
#      self.check_vtk("oxley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_oxley_3D_ReducedFunction(self):
#      dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ReducedFunction(dom).getX()
#      self.check_vtk("oxley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_oxley_3D_FunctionOnBoundary(self):
#      dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=FunctionOnBoundary(dom).getX()
#      self.check_vtk("oxley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

#   def test_oxley_3D_ReducedFunctionOnBoundary(self):
#      dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
#      x=ReducedFunctionOnBoundary(dom).getX()
#      self.check_vtk("oxley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
#                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
    
