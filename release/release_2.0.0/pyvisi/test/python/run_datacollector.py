
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
__url__="https://launchpad.net/escript-finley"

from esys.pyvisi import DataCollector
from esys.pyvisi.constant import *
from esys.escript import getMPISizeWorld
import unittest, os, sys

try:
	PYVISI_WORKDIR=os.environ['PYVISI_WORKDIR']
except KeyError:
	PYVISI_WORKDIR='.'
try:
	PYVISI_TEST_DATA_ROOT=os.environ['PYVISI_TEST_DATA_ROOT']
except KeyError:
	PYVISI_TEST_DATA_ROOT='.'

PYVISI_TEST_MESHES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, "data_meshes")
FILE_2D = "interior_2D.xml"
FILE_3D = "interior_3D.xml"

SCALAR_FIELD_POINT_DATA = "temperature"
VECTOR_FIELD_POINT_DATA = "velocity"
TENSOR_FIELD_POINT_DATA = "stress"
		
SCALAR_FIELD_CELL_DATA = "temperature_cell"
VECTOR_FIELD_CELL_DATA = "velocity_cell"
TENSOR_FIELD_CELL_DATA = "stress_cell"

class TestSourceXml:
	def setFileName(self, file):
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, file))

	def setActiveScalar(self, s):
		self.data_collector.setActiveScalar(scalar = s)

	def setActiveVector(self, v):
		self.data_collector.setActiveVector(vector = v)

	def setActiveTensor(self, t):
		self.data_collector.setActiveTensor(tensor = t)
	
	def checkScalarFieldPointData(self, f):
		self.failUnlessEqual(self.data_collector._getDataCollectorOutput().GetPointData().GetScalars().GetName(), f)
		
	def checkVectorFieldPointData(self, f):
		self.failUnlessEqual(self.data_collector._getDataCollectorOutput().GetPointData().GetVectors().GetName(), f)

	def checkTensorFieldPointData(self, f):
		self.failUnlessEqual(self.data_collector._getDataCollectorOutput().GetPointData().GetTensors().GetName(), f)

	def checkScalarFieldCellData(self, f):
		self.failUnlessEqual(self.data_collector._getDataCollectorOutput().GetCellData().GetScalars().GetName(), f)
		
	def checkVectorFieldCellData(self, f):
		self.failUnlessEqual(self.data_collector._getDataCollectorOutput().GetCellData().GetVectors().GetName(), f)

	def checkTensorFieldCellData(self, f):
		self.failUnlessEqual(self.data_collector._getDataCollectorOutput().GetCellData().GetTensors().GetName(), f)


class TestSourceXml2DPointData(unittest.TestCase, TestSourceXml):
	def setUp(self):
		self.data_collector = DataCollector(source = Source.XML)
		self.setFileName(FILE_2D)

	def tearDown(self):
		del self.data_collector

	def testSetActiveScalar(self):
		self.setActiveScalar(SCALAR_FIELD_POINT_DATA)
		self.checkScalarFieldPointData(SCALAR_FIELD_POINT_DATA)

	def testSetActiveVector(self):
		self.setActiveVector(VECTOR_FIELD_POINT_DATA)
		self.checkVectorFieldPointData(VECTOR_FIELD_POINT_DATA)

	def testSetActiveTensor(self):
		self.setActiveTensor(TENSOR_FIELD_POINT_DATA)
		self.checkTensorFieldPointData(TENSOR_FIELD_POINT_DATA)

class TestSourceXml2DCellData(unittest.TestCase, TestSourceXml):
	def setUp(self):
		self.data_collector = DataCollector(source = Source.XML)
		self.setFileName(FILE_2D)

	def tearDown(self):
		del self.data_collector

	def testSetActiveScalar(self):
		self.setActiveScalar(SCALAR_FIELD_CELL_DATA)
		self.checkScalarFieldCellData(SCALAR_FIELD_CELL_DATA)

	def testSetActiveVector(self):
		self.setActiveVector(VECTOR_FIELD_CELL_DATA)
		self.checkVectorFieldCellData(VECTOR_FIELD_CELL_DATA)

	def testSetActiveTensor(self):
		self.setActiveTensor(TENSOR_FIELD_CELL_DATA)
		self.checkTensorFieldCellData(TENSOR_FIELD_CELL_DATA)

class TestSourceXml3DPointData(unittest.TestCase, TestSourceXml):
	def setUp(self):
		self.data_collector = DataCollector(source = Source.XML)
		self.setFileName(FILE_3D)

	def tearDown(self):
		del self.data_collector

	def testSetActiveScalar(self):
		self.setActiveScalar(SCALAR_FIELD_POINT_DATA)
		self.checkScalarFieldPointData(SCALAR_FIELD_POINT_DATA)

	def testSetActiveVector(self):
		self.setActiveVector(VECTOR_FIELD_POINT_DATA)
		self.checkVectorFieldPointData(VECTOR_FIELD_POINT_DATA)

	def testSetActiveTensor(self):
		self.setActiveTensor(TENSOR_FIELD_POINT_DATA)
		self.checkTensorFieldPointData(TENSOR_FIELD_POINT_DATA)

class TestSourceXml3DCellData(unittest.TestCase, TestSourceXml):
	def setUp(self):
		self.data_collector = DataCollector(source = Source.XML)
		self.setFileName(FILE_3D)

	def tearDown(self):
		del self.data_collector

	def testSetActiveScalar(self):
		self.setActiveScalar(SCALAR_FIELD_CELL_DATA)
		self.checkScalarFieldCellData(SCALAR_FIELD_CELL_DATA)

	def testSetActiveVector(self):
		self.setActiveVector(VECTOR_FIELD_CELL_DATA)
		self.checkVectorFieldCellData(VECTOR_FIELD_CELL_DATA)

	def testSetActiveTensor(self):
		self.setActiveTensor(TENSOR_FIELD_CELL_DATA)
		self.checkTensorFieldCellData(TENSOR_FIELD_CELL_DATA)


###############################################################################


if __name__ == '__main__':
    if getMPISizeWorld() == 1: 
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml2DPointData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml2DCellData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml3DPointData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml3DCellData))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
    else:
        print "run_datacollector.py is not executed as more than one processor is used."

