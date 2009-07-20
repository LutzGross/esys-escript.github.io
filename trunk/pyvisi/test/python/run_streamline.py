
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
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

from esys.pyvisi import DataCollector, Scene, StreamLine, GlobalPosition
from esys.pyvisi.constant import *
import unittest, os, sys
from stat import ST_SIZE
from esys.escript import getMPISizeWorld

try:
	PYVISI_WORKDIR=os.environ['PYVISI_WORKDIR']
except KeyError:
	PYVISI_WORKDIR='.'
try:
	PYVISI_TEST_DATA_ROOT=os.environ['PYVISI_TEST_DATA_ROOT']
except KeyError:
	PYVISI_TEST_DATA_ROOT='.'

PYVISI_TEST_MESHES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, "data_meshes")
PYVISI_TEST_STREAMLINE_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_sample_images", "streamline")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestStreamLine:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_STREAMLINE_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(\
				PYVISI_TEST_STREAMLINE_IMAGES_PATH, file))[ST_SIZE] > \
				MIN_IMAGE_SIZE)

class TestStreamLinePointSource(unittest.TestCase, TestStreamLine):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.streamline = StreamLine(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.VECTOR,
				lut = Lut.COLOR, cell_to_point = False, outline = True)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.streamline

	def testPointSource(self):
		self.streamline.setPointSourceRadius(0.1)
		self.streamline.setPointSourceCenter(GlobalPosition(1.3, 1.3, 0.4))
		self.streamline.setPointSourceNumberOfPoints(2)
		self.render("TestStreamLinePointSource.jpg")

class TestStreamLineModule(unittest.TestCase, TestStreamLine):	
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.streamline = StreamLine(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.VECTOR,
				lut = Lut.COLOR, cell_to_point = False, outline = True)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.streamline

	def testStreamLineModule(self):
		self.streamline.setMaximumPropagationTime(20)
		self.streamline.setStepLength(0.1)
		self.streamline.setIntegrationStepLength(0.1)
		self.streamline.setIntegrationToBothDirections()
		self.render("TestStreamLineModule.jpg")

class TestStreamLineTube(unittest.TestCase, TestStreamLine):	
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.streamline = StreamLine(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.SCALAR,
				lut = Lut.COLOR, cell_to_point = False, outline = True)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.streamline

	def testSetTubeRadiusToVaryByVector(self):
		self.streamline.setTubeRadius(0.02)
		self.streamline.setTubeNumberOfSides(3)
		self.streamline.setTubeRadiusToVaryByVector()
		self.render("TestStreamLineTube_testSetTubeRadiusToVaryByVector.jpg")

	def testSetTubeRadiusToVaryByScalar(self):
		self.streamline.setTubeRadiusToVaryByScalar()
		self.render("TestStreamLineTube_testSetTubeRadiusToVaryByScalar.jpg")
		

##############################################################################


if __name__ == '__main__':
    if getMPISizeWorld() == 1: 
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLinePointSource))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineModule))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineTube))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
    else:
        print "run_streamline.py is not executed as more than one processor is used."

