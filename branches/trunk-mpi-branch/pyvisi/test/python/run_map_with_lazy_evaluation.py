#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

# Import the necessary modules
from esys.pyvisi import DataCollector, Scene, Map, MapOnPlaneCut, MapOnPlaneClip
from esys.pyvisi import Camera, MapOnScalarClip, MapOnScalarClipWithRotation
from esys.pyvisi.constant import *
import unittest, os
from stat import ST_SIZE

try:
	PYVISI_WORKDIR=os.environ['PYVISI_WORKDIR']
except KeyError:
	PYVISI_WORKDIR='.'
try:
	PYVISI_TEST_DATA_ROOT=os.environ['PYVISI_TEST_DATA_ROOT']
except KeyError:
	PYVISI_TEST_DATA_ROOT='.'

PYVISI_TEST_MESHES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, "data_meshes")
PYVISI_TEST_MAP_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "map")

MIN_IMAGE_SIZE = 100
FILE_3D_1 = "results.xml"
FILE_3D_2 = "interior_3D.xml"
FILE_2D_ROTATION = "without_st.0700.xml"
X_SIZE = 400
Y_SIZE = 400

JPG_RENDERER = Renderer.OFFLINE_JPG
SCALAR_FIELD_POINT_DATA_1 = "scalar1"
SCALAR_FIELD_POINT_DATA_2 = "scalar2"

class TestMapWithLazyEvaluation:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_MAP_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_MAP_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestMapLazy(unittest.TestCase, TestMapWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testMapLazy(self):

		# Create a scene with four viewports.
		s = Scene(renderer = JPG_RENDERER, num_viewport = 4, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)
		dc1.setActiveScalar(scalar = SCALAR_FIELD_POINT_DATA_1)

		dc2 = DataCollector(source = Source.XML)
		dc2.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D_1))
		dc2.setActiveScalar(scalar = SCALAR_FIELD_POINT_DATA_2)

		# Create a map instance for the first viewport.
		m1 = Map(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)

		# Create a map instance for the second viewport.
		m2 = Map(scene = s, data_collector = dc2, 
				viewport = Viewport.NORTH_WEST, lut = Lut.COLOR, outline = True)
		m2.setColor(color = Color.BLUE)

		# Create a map instance for the third viewport.
		m3 = Map(scene = s, data_collector = dc1, 
				viewport = Viewport.NORTH_EAST, lut = Lut.COLOR, outline = True)

		# Create a map instance the fourth viewport.
		m4 = Map(scene = s, data_collector = dc2, 
				viewport = Viewport.SOUTH_EAST, lut = Lut.COLOR, outline = True)

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D_1))
		m1.setRepresentationToWireframe()
		m4.setOpacity(opacity = 0.5)
	
		self.render("TestMapWithLazyEvaluation.jpg")

class TestMapOnPlaneCutLazy(unittest.TestCase, TestMapWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testMapOnPlaneCutLazy(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 4, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		# Create two data collector instances for two different sources.
		dc1 = DataCollector(source = Source.XML)

		# Create a map on plane cut instance for the first viewport.
		mopc1 = MapOnPlaneCut(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST)
		mopc1.setPlaneToYZ(offset = 1.5)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.isometricView()

		# Create three map on plane cut instances for the second viewport.
		mopc2_1 = MapOnPlaneCut(scene = s, data_collector = dc1, 
				viewport = Viewport.NORTH_WEST)

		mopc2_2 = MapOnPlaneCut(scene = s, data_collector = dc1, 
				viewport = Viewport.NORTH_WEST)

		mopc2_3 = MapOnPlaneCut(scene = s, data_collector = dc1, 
				viewport = Viewport.NORTH_WEST)
		mopc2_3.setPlaneToXY()
		mopc2_3.translate(0,0,0.5)

		c2 = Camera(scene = s, viewport = Viewport.NORTH_WEST)
		c2.isometricView()

		mopc2_1.setPlaneToYZ(offset = 1.5)
		mopc2_2.setPlaneToXZ(offset = 1.5)
		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D_2)) 

		self.render("TestMapOnPlaneCutWithLazyEvaluation.jpg")


class TestMapOnPlaneClipLazy(unittest.TestCase, TestMapWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testMapOnPlaneClipLazy(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create three map on clip instances.
		mopc1_1 = MapOnPlaneClip(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST)
		mopc1_1.setPlaneToXY()
		mopc1_1.rotateX(angle = 5)

		mopc1_2 = MapOnPlaneClip(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST)

		mopc1_3 = MapOnPlaneClip(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST)
		mopc1_3.setPlaneToXZ()
		mopc1_3.rotateX(angle = -40)
		mopc1_3.translate(x_offset = 0, y_offset = 2.2, z_offset = 0)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.isometricView()

		mopc1_2.setPlaneToYZ(offset = 2.5)
		mopc1_2.setOpacity(opacity = 0.5)
		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D_2))

		self.render("TestMapOnPlaneClipWithLazyEvaluation.jpg")


class TestMapOnScalarClipLazy(unittest.TestCase, TestMapWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testMapOnScalarClipLazy(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create a map on scalar clip instance.
		mosc1_1 = MapOnScalarClip(scene = s, data_collector = dc1, 
				lut = Lut.GREY_SCALE)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.isometricView()

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D_2))
		self.render("TestMapOnScalarClipWithLazyEvaluation.jpg")

class TestMapOnScalarClipWithRotationLazy(unittest.TestCase, 
		TestMapWithLazyEvaluation):

	def tearDown(self):
		del self.scene

	def testMapOnScalarClipWithRotationLazy(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create a map on scalar clip instance.
		mosc1_1 = MapOnScalarClipWithRotation(scene = s, data_collector = dc1, 
				lut = Lut.GREY_SCALE)
		mosc1_1.setAngle(200)
		
		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.isometricView()

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_2D_ROTATION))
		self.render("TestMapOnScalarClipWithRotationWithLazyEvaluation.jpg")





##############################################################################
if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneCutLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneClipLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClipLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClipWithRotationLazy))
	unittest.TextTestRunner(verbosity=2).run(suite)
