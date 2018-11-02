
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for input and output of meshes and data objects

:remark:

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

import esys.escriptcore.utestselect as unittest, sys
from esys.escriptcore.testing import *
from esys.escript import *
from esys.finley import Rectangle, Brick, LoadMesh, ReadMesh, GetMeshFromFile, ReadGmsh
import os as os

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR='.'

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA=os.path.join(os.getcwd(), "finley/test/python")

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

REL_TOL=1.e-6

# Number of elements scales up with number of MPI processes
mpisize = getMPISizeWorld()
NE0 = 7 * mpisize
NE1 = 11
NE2 = 5

class Test_InputOutputOnFinley(unittest.TestCase):
     # Check that two domains are equal using Fourier integrals
     # We cannot compare the X coordinates since they are on different domains
     def domainsEqual(self, m1, m2, nft=100):
        self.assertTrue(m1.getDim() == m2.getDim(), "Dimensions differ")
        self.assertTrue(m1.getNumDataPointsGlobal() == m2.getNumDataPointsGlobal(), "Num data points differ")
        for tagName in m1.showTagNames().split(", "):
          self.assertTrue(m2.isValidTagName(tagName), "m1 has a tag '%s' not present in m2" % tagName)
        for tagName in m2.showTagNames().split(", "):
          self.assertTrue(m1.isValidTagName(tagName), "m2 has a tag '%s' not present in m1" % tagName)
          self.assertTrue(m1.getTag(tagName) == m2.getTag(tagName), "values of tag '%s' differ" % tagName)
        for fs in ["Solution", "ReducedSolution", "Function", "ReducedFunction", "ContinuousFunction", "ReducedContinuousFunction"]:
          fs1 = eval("%s(m1)" % fs)
          fs2 = eval("%s(m2)" % fs)
          x1 = fs1.getX()
          x2 = fs2.getX()
          for n in range(1, nft+1):
            integ1 = integrate(sin(n*x1))
            integ2 = integrate(sin(n*x2))
            self.assertTrue(Lsup(abs(integ1-integ2)) <= REL_TOL, "integrals for n=%d differ" % n)
        return True

     # Does optimize=True change Rectangle for order=1?
     def test_Rectangle_optimize_order1(self):
        mydomain1 = Rectangle(n0=NE0, n1=NE1, order=1, l0=1., l1=1., optimize=False,useElementsOnFace=0)
        mydomain2 = Rectangle(n0=NE0, n1=NE1, order=1, l0=1., l1=1., optimize=True,useElementsOnFace=0)
        self.domainsEqual(mydomain1, mydomain2)

     # Does optimize=True change Rectangle for order=2?
     def test_Rectangle_optimize_order2(self):
        mydomain1 = Rectangle(n0=NE0, n1=NE1, order=2, l0=1., l1=1., optimize=False,useElementsOnFace=0)
        mydomain2 = Rectangle(n0=NE0, n1=NE1, order=2, l0=1., l1=1., optimize=True,useElementsOnFace=0)
        self.domainsEqual(mydomain1, mydomain2)

     # Does optimize=True change Rectangle for order=-1?
     def test_Rectangle_optimize_macro(self):
        mydomain1 = Rectangle(n0=NE0, n1=NE1, order=-1, l0=1., l1=1., optimize=False, useElementsOnFace=0)
        mydomain2 = Rectangle(n0=NE0, n1=NE1, order=-1, l0=1., l1=1., optimize=True,useElementsOnFace=0)
        self.domainsEqual(mydomain1, mydomain2)

     # Does optimize=True change Brick for order=1?
     def test_Brick_optimize_order1(self):
        mydomain1 = Brick(n0=NE0, n1=NE1, n2=NE2, order=1, l0=1., l1=1., l2=1., optimize=False,useElementsOnFace=0)
        mydomain2 = Brick(n0=NE0, n1=NE1, n2=NE2, order=1, l0=1., l1=1., l2=1., optimize=True,useElementsOnFace=0)
        self.domainsEqual(mydomain1, mydomain2)

     # Does optimize=True change Brick for order=2?
     def test_Brick_optimize_order2(self):
        mydomain1 = Brick(n0=NE0, n1=NE1, n2=NE2, order=2, l0=1., l1=1., l2=1., optimize=False,useElementsOnFace=0)
        mydomain2 = Brick(n0=NE0, n1=NE1, n2=NE2, order=2, l0=1., l1=1., l2=1., optimize=True,useElementsOnFace=0)
        self.domainsEqual(mydomain1, mydomain2)
     # Does optimize=True change Brick for order=-1?
     def test_Brick_optimize_macro(self):
        mydomain1 = Brick(n0=NE0, n1=NE1, n2=NE2, order=-1, l0=1., l1=1., l2=1., optimize=False,useElementsOnFace=0)
        mydomain2 = Brick(n0=NE0, n1=NE1, n2=NE2, order=-1, l0=1., l1=1., l2=1., optimize=True,useElementsOnFace=0)
        self.domainsEqual(mydomain1, mydomain2)

     @unittest.skipIf(not loadIsConfigured(), "load not configured")
     def test_data_dump_to_NetCDF_rectangle(self):
        mydomain1 = Rectangle(n0=NE0, n1=NE1, order=1, l0=1., l1=1., optimize=False,useElementsOnFace=0)
        d1=Data(mydomain1.getMPIRank(), Function(mydomain1))
        d1.expand()
        dumpfile=os.path.join(FINLEY_WORKDIR, "tempfile.dump.nc")
        d1.dump(dumpfile)
        d2=load(dumpfile, mydomain1)
        self.assertTrue(Lsup(abs(d1-d2)) <= REL_TOL, "data objects differ")

     @unittest.skipIf(not loadIsConfigured(), "load not configured")
     def test_data_dump_to_NetCDF_brick(self):
        mydomain1 = Brick(n0=NE0, n1=NE1, n2=NE2, order=2, l0=1., l1=1., l2=1., optimize=False)
        d1=Data(mydomain1.getMPIRank(), Function(mydomain1))
        d1.expand()
        dumpfile=os.path.join(FINLEY_WORKDIR, "tempfile.dump.nc")
        d1.dump(dumpfile)
        d2=load(dumpfile, mydomain1)
        self.assertTrue(Lsup(abs(d1-d2)) <= REL_TOL, "data objects differ")

     @unittest.skipIf(not loadIsConfigured() or mpisize > 1, "load not configured or mpisize > 1")
     def test_mesh_dump_to_NetCDF_rectangle(self):
        mydomain1 = Rectangle(n0=NE0, n1=NE1, order=1, l0=1., l1=1., optimize=False)
        dumpfile=os.path.join(FINLEY_WORKDIR, "tempfile.mesh.nc")
        mydomain1.dump(dumpfile)
        mydomain2=LoadMesh(dumpfile)
        self.domainsEqual(mydomain1, mydomain2)

     @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
     def test_gmshTags(self):
        dom=ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH, "tagtest.msh"),2)
        tags=dom.showTagNames().split(', ')
        self.assertEqual(tags,['tag1', 'tag2', 'tag3'],'error with tags')
        self.assertEqual(dom.getTag('tag1'),1,'error with tag1')
        self.assertEqual(dom.getTag('tag2'),2,'error with tag2')
        self.assertEqual(dom.getTag('tag3'),3,'error with tag3')
        self.assertRaises(ValueError, dom.getTag, 'tag4')
     
     @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
     def test_gmshNamedTags(self):
        dom=ReadGmsh(os.path.join(FINLEY_TEST_MESH_PATH, "test_Add.msh"),3)
        sigma = Scalar(0,Function(dom))
        sigma.expand()
        sigma.setTaggedValue("v5",10)
        self.assertEqual(Lsup(sigma),10.0,"Named value not set properly")

     def test_flyTags(self):
        dom=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH, "rectangle_8x10.fly"))
        tags=sorted(dom.showTagNames().split(', '))
        self.assertEqual(tags,sorted(['top', 'bottom', 'left', 'right']))
        self.assertEqual(dom.getTag('top'),20,'error with top')
        self.assertEqual(dom.getTag('bottom'),10,'error with bottom,')
        self.assertEqual(dom.getTag('left'),1,'error with left')
        self.assertEqual(dom.getTag('right'),2,'error with reight')
        self.assertRaises(ValueError, dom.getTag, 'tag4')

     @unittest.skipIf(not loadIsConfigured() or mpisize>1, "load not configured or mpisize > 1")
     def test_mesh_dump_to_NetCDF_brick(self):
        mydomain1 = Brick(n0=NE0, n1=NE1, n2=NE2, order=2, l0=1., l1=1., l2=1., optimize=False)
        dumpfile=os.path.join(FINLEY_WORKDIR, "tempfile.mesh.nc")
        mydomain1.dump(dumpfile)
        mydomain2=LoadMesh(dumpfile)
        self.domainsEqual(mydomain1, mydomain2)

     @unittest.skipIf(mpisize>15, "more than 15 MPI ranks")
     def test_mesh_read_rectangle_from_finley_file(self):
         mydomain1 = Rectangle(n0=8, n1=10, order=1, l0=1., l1=1., optimize=False)
         mydomain2 = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"rectangle_8x10.fly"))
         self.domainsEqual(mydomain1, mydomain2)

     @unittest.skipIf(mpisize>15, "more than 15 MPI ranks")
     def test_mesh_read_brick_from_finley_file(self):
          mydomain1 = Brick(n0=8, n1=10, n2=12, order=1, l0=1., l1=1., l2=1., optimize=False)
          mydomain2 = ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"brick_8x10x12.fly"))
          self.domainsEqual(mydomain1, mydomain2)

     @unittest.skipIf(mpisize>1, "more than 1 MPI rank")
     def test_GetMeshFromFile(self):
           m=GetMeshFromFile(os.path.join(FINLEY_TEST_MESH_PATH,'tet10_gmsh.msh'), numDim=3)
           del m
           m=GetMeshFromFile(os.path.join(FINLEY_TEST_MESH_PATH, 'tet10.fly'))
           # now we try some params
           m=GetMeshFromFile(os.path.join(FINLEY_TEST_MESH_PATH,'tet10_gmsh.msh'), numDim=3, integrationOrder=2)

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

