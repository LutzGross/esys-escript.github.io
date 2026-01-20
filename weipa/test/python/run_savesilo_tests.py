
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
from esys.escript import ContinuousFunction, Function, ReducedFunction,\
            FunctionOnBoundary, ReducedFunctionOnBoundary,\
            FunctionOnContactZero, ReducedFunctionOnContactZero,\
            FunctionOnContactOne, ReducedFunctionOnContactOne,\
            Solution, ReducedSolution, getMPISizeWorld
from esys.weipa import saveSilo

try:
    import Silo
    HAVE_SILO=True
except ImportError:
    HAVE_SILO=False
try:
    from esys import finley
    HAVE_FINLEY=True
except ImportError:
    HAVE_FINLEY=False
try:
    from esys import ripley
    HAVE_RIPLEY=True
except ImportError:
    HAVE_RIPLEY=False
try:
    from esys import oxley
    HAVE_OXLEY=True
except ImportError:
    HAVE_OXLEY=False

try:
    WEIPA_TEST_MESHES=os.environ['WEIPA_TEST_MESHES']
except KeyError:
    WEIPA_TEST_MESHES='meshes'

try:
    WEIPA_WORKDIR=os.environ['WEIPA_WORKDIR']
except KeyError:
    WEIPA_WORKDIR=os.path.join(os.getcwd(),'weipa/test/python/')

class SiloReader():
    """
    Silo file reader that uses the Python interface to the Silo library.
    """
    def __init__(self):
        self.f=None

    def open(self, filename):
        try:
            self.f=Silo.Open(filename)
            self.filename=filename
            return True
        except:
            return False

    def close(self):
        if self.f:
            self.f.Close()

    def getDirNames(self):
        return self.f.GetToc().dir_names

    def setDir(self, dirname):
        self.f.SetDir(dirname)

    def getTimeAndCycle(self):
        return self.f.GetVar('dtime'), self.f.GetVar('cycle')

    def getElementNames(self):
        names=[]
        for v in self.f.GetToc().var_names:
            if v.endswith('_coord0'):
                names.append(v[:-7])
        return names

    def getCellCoords(self, elementName):
        ret=()
        for i in range(3):
            try:
                data=self.f.GetVar("%s_coord%d"%(elementName,i))
                ret+=data,
            except:
                ret+=None,
        return ret

    def getCellExtents(self, elementName):
        ret=()
        for i in ('min','max'):
            try:
                data=self.f.GetVar("%s_%s_extents"%(elementName,i))
                ret+=data,
            except:
                ret+=None,
        return ret

    def getCellNodeList(self, elementName):
        try:
            return self.f.GetVar("%s_zones_nodelist"%elementName)
        except:
            return None

    def getCellShapeInfo(self, elementName):
        infoNames=('zones_shapecnt', 'zones_shapesize', 'zones_shapetype')
        ret=()
        for n in infoNames:
            try:
                data=self.f.GetVar("%s_%s"%(elementName,n))
                ret+=data,
            except:
                ret+=None,
        return ret

    def getData(self):
        data={}
        names=self.f.GetToc().var_names
        for v in names:
            if v[-5:]=='_data':
                data[v[:-5]] = self.f.GetVar(v)
        return data

class SiloSaver(unittest.TestCase): #requires subclassing
    def numericCompareL2(self, vector1, vector2):
        """
        Compares two lists of numbers using the L2 norm, returns true if they
        match up to a tolerance TOL, false otherwise
        """
        if vector2 == None: return False
        TOL = 2.0e-5
        try:
            l1=len(vector1)
        except TypeError:
            vector1=[vector1]
        try:
            l2=len(vector2)
        except TypeError:
            vector2=[vector2]
        if len(vector1) != len(vector2): return False
        diff = 0.0
        for i in range(0, len(vector1)):
            tmp = vector1[i] - vector2[i]
            diff += tmp * tmp
        if math.sqrt(diff) > TOL: return False
        return True

    def compareSiloFiles(self, file1, file2):
        """
        Compares two Silo files, asserts if any element is not equal.
        file2 is the reference file to compare against
        """
        p1=SiloReader()
        p2=SiloReader()
        self.assertTrue(p1.open(file1), "Invalid Silo file '%s'"%file1)
        p2.open(file2)
        self.assertEqual(p1.getTimeAndCycle(), p2.getTimeAndCycle())

        dirNames1=p1.getDirNames()
        dirNames2=p2.getDirNames()
        for dirName in dirNames2:
            self.assertTrue(dirName in dirNames1, "Silo file is missing directories")
        p1.setDir('/block0000')
        p2.setDir('/block0000')

        elementNames1=p1.getElementNames()
        elementNames2=p2.getElementNames()

        for elementName in elementNames2:
            self.assertTrue(elementName in elementNames1, "Mesh '%s' missing in Silo file"%elementName)

            cx1,cy1,cz1=p1.getCellCoords(elementName)
            cx2,cy2,cz2=p2.getCellCoords(elementName)
            self.assertEqual(len(cx1), len(cx2))
            self.assertEqual(len(cy1), len(cy2))
            if cz2 is not None:
                self.assertEqual(len(cz1), len(cz2))
                coords1=zip(cx1,cy1,cz1)
                coords2=zip(cx2,cy2,cz2)
            else:
                coords1=zip(cx1,cy1)
                coords2=zip(cx2,cy2)

            # Find mapping of nodes in file 1 to file 2 (they may be
            # permuted)
            #nodeMap1to2 = {}
            for i1 in range(len(coords1)):
                indexList=[]
                for i2 in range(len(coords2)):
                    if self.numericCompareL2(coords1[i1], coords2[i2]):
                        indexList.append(i2)
                self.assertNotEquals(len(indexList), 0,
                                 "Node with coordinates %s missing in '%s'"%(str(coords1[i1]),elementName))
                #nodeMap1to2[i1]=indexList

            extents1=p1.getCellExtents(elementName)
            extents2=p2.getCellExtents(elementName)
            self.assertEqual(extents1, extents2)

            nodelist1=p1.getCellNodeList(elementName)
            nodelist2=p2.getCellNodeList(elementName)

            ccount1,csize1,ctype1=p1.getCellShapeInfo(elementName)
            ccount2,csize2,ctype2=p2.getCellShapeInfo(elementName)

        # data
        data1=p1.getData()
        data2=p2.getData()
        p1.close()
        p2.close()

        # TODO: The Silo module does not allow checking for the mesh of
        # a variable yet so we cannot compare permuted entries using the
        # node mappings from above (see vtk tests)
        self.assertEqual(len(data1), len(data2))
        for name in data2:
            self.assertTrue(name in data1, "Variable '%s' missing"%name)
            if name.startswith('data'):
                self.assertTrue(self.numericCompareL2(
                    data1[name], data2[name]),
                    "Values in '%s' do not match" % name)

    def check_silo(self, reference, **data):
        outFileBase="out_"+reference
        saveSilo(os.path.join(WEIPA_WORKDIR, outFileBase), write_meshdata=True, **data)
        ref=os.path.join(WEIPA_TEST_MESHES, reference+".silo")
        out=os.path.join(WEIPA_WORKDIR, outFileBase+".silo")
        self.compareSiloFiles(out, ref)
        self.cleanup(out)

    def cleanup(self, filename):
        os.remove(filename)

class Test_Silo_import(unittest.TestCase):
    def test_import(self):
        if not HAVE_SILO:
            try:
                import Silo
            except ImportError as e:
                if "No module named Silo" not in str(e):
                    raise unittest.SkipTest("Silo module broken")

@unittest.skipIf(getMPISizeWorld()>1, "MPI size > 1")
@unittest.skipIf(not HAVE_FINLEY, "finley module not available")
@unittest.skipIf(not HAVE_SILO, "Silo module not available")
class Test_Finley_SaveSilo(SiloSaver):

  # === Finley hex 2D order 1 with contacts ===================================

  def test_hex_contact_2D_order1_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_onFace_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                             data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_onFace_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                             data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_onFace_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_onFace_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_onFace_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order1_onFace_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order1_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Finley hex 2D order 2 with contacts ===================================

  def test_hex_contact_2D_order2_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_2D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_2D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_2D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_onFace_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                             data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_onFace_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                             data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o2_rcontact", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_onFace_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                            data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_onFace_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                            data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o2_rcontact", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_onFace_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                            data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_contact_2D_order2_onFace_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_2D_order2_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_2D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.],
                                            data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Finley hex 2D order 2 =================================================

  def test_hex_2D_order2_empty(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
     self.check_silo("hex_2D_o2", domain=dom)

  def test_hex_2D_order2_AllPoints_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o2_node_3xs", data_r=x_r[0], data_n=x_n[0], data=x[0])

  def test_hex_2D_order2_2Cells_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2_cell_2xs", data_b=x_b[0], data=x[0])

  def test_hex_2D_order2_BoundaryPoint_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2_boundary_2xs", data=x[0],data_b=x_b[0])

  def test_hex_2D_order2_Cells_AllData(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_2D_o2_cell_all",
                     data_s=x[0],
                     data_v=x[0]*[1.,2.],
                     data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2_CellsPoints_AllData(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2.msh"),optimize=False)
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o2_cellnode_all",
                     data_sp=x_p[0],
                     data_vp=x_p[0]*[1.,2.],
                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                     data_sc=x_c[0],
                     data_vc=x_c[0]*[1.,2.],
                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

  # === Finley hex 2D order 2 (full) ==========================================

  def test_hex_2D_order2p_empty(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     self.check_silo("hex_2D_o2p", domain=dom)

  def test_hex_2D_order2p_AllPoints_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o1_node_3xs", data_r=x_r[0], data_n=x_n[0], data=x[0])

  def test_hex_2D_order2p_2Cells_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2p_cell_2xs", data_b=x_b[0], data=x[0])

  def test_hex_2D_order2p_BoundaryPoint_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2p_boundary_2xs", data=x[0],data_b=x_b[0])

  def test_hex_2D_order2p_CellsPoints_AllData(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o2p_cellnode_all",
                     data_sp=x_p[0],
                     data_vp=x_p[0]*[1.,2.],
                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                     data_sc=x_c[0],
                     data_vc=x_c[0]*[1.,2.],
                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2p_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2p_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2p_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_2D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
                                         data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2p_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_2D_o2p_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2p_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_2D_o2p_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
                                         data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2p_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2p_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                            data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_order2p_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_order2p.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_o2p_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                            data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Finley hex 2D macro ===================================================

  def test_hex_2D_macro_empty(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     self.check_silo("hex_2D_o2p", domain=dom)

  def test_hex_2D_macro_AllPoints_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o1_node_3xs", data_r=x_r[0], data_n=x_n[0], data=x[0])

  def test_hex_2D_macro_CellsPoints(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_macro_cellnode_all",
                     data_sp=x_p[0],
                     data_vp=x_p[0]*[1.,2.],
                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                     data_sc=x_c[0],
                     data_vc=x_c[0]*[1.,2.],
                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_macro_2Cells_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_cell_2xs", data_b=x_b[0], data=x[0])

  def test_hex_2D_macro_BoundaryPoint_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_boundary_2xs", data_b=x_b[0], data=x[0])

  def test_hex_2D_macro_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_macro_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_2D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_macro_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_2D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
                                         data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_macro_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_2D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_macro_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_2D_macro_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_macro_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                              data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_hex_2D_macro_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_2D_macro.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                              data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Finley hex 3D order 1 with contacts ===================================

  def test_hex_contact_3D_order1_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_onFace_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_onFace_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o1_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_onFace_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_onFace_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_onFace_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order1_onFace_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order1_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o1_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  # === Finley hex 3D order 2 with contacts ===================================

  def test_hex_contact_3D_order2_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_3D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_3D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_3D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_onFace_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_onFace_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o2_f_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_onFace_FunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
     x=FunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_onFace_ReducedFunctionOnContactZero(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactZero(dom).getX()
     self.check_silo("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o2_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_onFace_FunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
     x=FunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_contact_3D_order2_onFace_ReducedFunctionOnContactOne(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_contact_3D_order2_onFace.msh"),optimize=False)
     x=ReducedFunctionOnContactOne(dom).getX()
     self.check_silo("hex_3D_o2_f_contact", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  # === Finley hex 3D order 2 (full) ==========================================

  def test_hex_3D_order2p_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_order2p_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_order2p_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_3D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                         data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_order2p_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_3D_o2p_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_order2p_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_3D_o2p_rcell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                         data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_order2p_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o2p_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_order2p_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_order2p.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_o2p_rboundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                             data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  # === Finley hex 3D macro ===================================================

  def test_hex_3D_macro_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_macro_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("hex_3D_o2p_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_macro_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("hex_3D_o2p_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                         data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_macro_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("hex_3D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_macro_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("hex_3D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_macro_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                              data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_hex_3D_macro_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"hex_3D_macro.msh"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("hex_3D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                              data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  # === Finley tet 2D order 1 =================================================

  def test_tet_2D_order1_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order1_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("tet_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order1_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("tet_2D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order1_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("tet_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order1_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("tet_2D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order1_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order1_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order1.fly"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Finley tet 2D order 2 =================================================

  def test_tet_2D_order2(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     self.check_silo("tet_2D_o2", domain=dom)

  def test_tet_2D_order2_AllPoints_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o1_node_3xs", data_r=x_r[0], data_n=x_n[0], data=x[0])

  def test_tet_2D_order2_02Points_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=Solution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o2_node_2xs", data_n=x_n[0], data=x[0])

  def test_tet_2D_order2_2Cells_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_o2_cell_2xs", data_b=x_b[0], data=x[0])

  def test_tet_2D_order2_BoundaryPoint_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_o2_boundary_2xs", data=x[0],data_b=x_b[0])

  def test_tet_2D_order2_Cells_AllData(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("tet_2D_o2_cell_all",
                     data_s=x[0],
                     data_v=x[0]*[1.,2.],
                     data_t=x[0]*[[11.,12.],[21.,22.]],
                     data_t2=x[0]*[[-11.,-12.],[-21.,-22.]])

  def test_tet_2D_order2_CellsPoints_AllData(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o2_cellnode_all",
                     data_sp=x_p[0],
                     data_vp=x_p[0]*[1.,2.],
                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                     data_sc=x_c[0],
                     data_vc=x_c[0]*[1.,2.],
                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order2_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order2_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order2_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("tet_2D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order2_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("tet_2D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order2_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("tet_2D_o2_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order2_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_order2_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_order2.fly"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Finley tet 2D macro ===================================================

  def test_tet_2D_macro(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     self.check_silo("tet_2D_o2", domain=dom)

  def test_tet_2D_macro_AllPoints_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=Solution(dom).getX()
     x_r=ReducedSolution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o1_node_3xs", data_r=x_r[0], data_n=x_n[0], data=x[0])

  def test_tet_2D_macro_02Points_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=Solution(dom).getX()
     x_n=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o2_node_2xs", data_n=x_n[0], data=x[0])

  def test_tet_2D_macro_2Cells_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=Function(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_cell_2xs", data_b=x_b[0], data=x[0])

  def test_tet_2D_macro_BoundaryPoint_Scalar(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     x_b=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_boundary_2xs", data_b=x_b[0], data=x[0])

  def test_tet_2D_macro_CellsPoints_AllData(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x_c=Function(dom).getX()
     x_p=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_macro_cellnode_all",
                     data_sp=x_p[0],
                     data_vp=x_p[0]*[1.,2.],
                     data_tp=x_p[0]*[[11.,12.],[21.,22.]],
                     data_sc=x_c[0],
                     data_vc=x_c[0]*[1.,2.],
                     data_tc=x_c[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_macro_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_macro_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("tet_2D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_macro_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("tet_2D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.],
                                        data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_macro_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("tet_2D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                          data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_macro_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("tet_2D_macro_rcell", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_macro_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                              data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_tet_2D_macro_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_2D_macro.fly"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("tet_2D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                              data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Finley tet 3D order 1 =================================================

  def test_tet_3D_order1_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("tet_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order1_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("tet_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order1_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("tet_3D_o1_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order1_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("tet_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order1_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("tet_3D_o1_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order1_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order1_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order1.fly"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("tet_3D_o1_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  # === Finley tet 3D order 2 =================================================

  def test_tet_3D_order2_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order2_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order2_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("tet_3D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order2_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("tet_3D_o2_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order2_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("tet_3D_o2_rcell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order2_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_3D_o2_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_order2_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_order2.fly"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("tet_3D_o2_rboundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                            data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  # === Finley tet 3D macro ===================================================

  def test_tet_3D_macro_ContinuousFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
     x=ContinuousFunction(dom).getX()
     self.check_silo("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_macro_Solution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
     x=Solution(dom).getX()
     self.check_silo("tet_3D_o2_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_macro_ReducedSolution(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
     x=ReducedSolution(dom).getX()
     self.check_silo("tet_3D_o2_rnode", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                        data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_macro_Function(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
     x=Function(dom).getX()
     self.check_silo("tet_3D_macro_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                          data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_macro_ReducedFunction(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
     x=ReducedFunction(dom).getX()
     self.check_silo("tet_3D_macro_rcell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_macro_FunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("tet_3D_macro_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                              data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_tet_3D_macro_ReducedFunctionOnBoundary(self):
     dom=finley.ReadMesh(os.path.join(WEIPA_TEST_MESHES,"tet_3D_macro.fly"),optimize=False)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("tet_3D_macro_rboundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                               data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

@unittest.skipIf(getMPISizeWorld()>1, "MPI size > 1")
@unittest.skipIf(not HAVE_SILO, "Silo module not available")
class Test_Ripley_SaveSilo(SiloSaver):

  # === Ripley 2D =============================================================

  def test_ripley_2D_ContinuousFunction(self):
     dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ContinuousFunction(dom).getX()
     self.check_silo("ripley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_ripley_2D_Solution(self):
     dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=Solution(dom).getX()
     self.check_silo("ripley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_ripley_2D_ReducedSolution(self):
     dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ReducedSolution(dom).getX()
     self.check_silo("ripley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_ripley_2D_Function(self):
     dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=Function(dom).getX()
     self.check_silo("ripley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_ripley_2D_ReducedFunction(self):
     dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ReducedFunction(dom).getX()
     self.check_silo("ripley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_ripley_2D_FunctionOnBoundary(self):
     dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("ripley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_ripley_2D_ReducedFunctionOnBoundary(self):
     dom=ripley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("ripley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  # === Ripley 3D =============================================================

  def test_ripley_3D_ContinuousFunction(self):
     dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ContinuousFunction(dom).getX()
     self.check_silo("ripley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_ripley_3D_Solution(self):
     dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=Solution(dom).getX()
     self.check_silo("ripley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_ripley_3D_ReducedSolution(self):
     dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ReducedSolution(dom).getX()
     self.check_silo("ripley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_ripley_3D_Function(self):
     dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=Function(dom).getX()
     self.check_silo("ripley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_ripley_3D_ReducedFunction(self):
     dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ReducedFunction(dom).getX()
     self.check_silo("ripley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_ripley_3D_FunctionOnBoundary(self):
     dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("ripley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_ripley_3D_ReducedFunctionOnBoundary(self):
     dom=ripley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("ripley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

@unittest.skipIf(getMPISizeWorld()>1, "MPI size > 1")
@unittest.skipIf(not HAVE_OXLEY, "ripley module not available")
@unittest.skipIf(not HAVE_SILO, "Silo module not available")
class Test_Oxley_SaveSilo(SiloSaver):

  # === Oxley 2D =============================================================

  def test_oxley_2D_ContinuousFunction(self):
     dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ContinuousFunction(dom).getX()
     self.check_silo("oxley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_oxley_2D_Solution(self):
     dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=Solution(dom).getX()
     self.check_silo("oxley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_oxley_2D_ReducedSolution(self):
     dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ReducedSolution(dom).getX()
     self.check_silo("oxley_2D_node", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_oxley_2D_Function(self):
     dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=Function(dom).getX()
     self.check_silo("oxley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_oxley_2D_ReducedFunction(self):
     dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ReducedFunction(dom).getX()
     self.check_silo("oxley_2D_cell", data_s=x[0], data_v=x[0]*[1.,2.],
                                       data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_oxley_2D_FunctionOnBoundary(self):
     dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("oxley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  def test_oxley_2D_ReducedFunctionOnBoundary(self):
     dom=oxley.Rectangle(n0=11, n1=3, l0=(-2.5,8.0), l1=(1.2,3.8), d0=getMPISizeWorld())
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("oxley_2D_boundary", data_s=x[0], data_v=x[0]*[1.,2.],
                                           data_t=x[0]*[[11.,12.],[21.,22.]])

  # === oxley 3D =============================================================

  def test_oxley_3D_ContinuousFunction(self):
     dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ContinuousFunction(dom).getX()
     self.check_silo("oxley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_oxley_3D_Solution(self):
     dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=Solution(dom).getX()
     self.check_silo("oxley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_oxley_3D_ReducedSolution(self):
     dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ReducedSolution(dom).getX()
     self.check_silo("oxley_3D_node", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_oxley_3D_Function(self):
     dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=Function(dom).getX()
     self.check_silo("oxley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_oxley_3D_ReducedFunction(self):
     dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ReducedFunction(dom).getX()
     self.check_silo("oxley_3D_cell", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                       data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_oxley_3D_FunctionOnBoundary(self):
     dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=FunctionOnBoundary(dom).getX()
     self.check_silo("oxley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])

  def test_oxley_3D_ReducedFunctionOnBoundary(self):
     dom=oxley.Brick(n0=11, n1=3, n2=2, l0=(-2.5,7.0), l1=(1.2,3.8), l2=4., d0=getMPISizeWorld(), d1=1, d2=1)
     x=ReducedFunctionOnBoundary(dom).getX()
     self.check_silo("oxley_3D_boundary", data_s=x[0], data_v=x[0]*[1.,2.,3.],
                                           data_t=x[0]*[[11.,12.,13.],[21.,22.,23.],[31.,32.,33.]])


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

