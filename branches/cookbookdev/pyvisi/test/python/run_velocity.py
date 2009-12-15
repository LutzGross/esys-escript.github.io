
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

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.pyvisi import DataCollector, Scene, Velocity
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
PYVISI_TEST_VELOCITY_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_sample_images", "velocity")

MIN_IMAGE_SIZE = 100
FILE_2D = "interior_2D.xml"
FILE_3D = "interior_3D.xml"
FILE_SECOND_ORDER_3D = "vel-000719.vtu"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestVelocity:
	def render(self, file):
		self.scene.render(image_name = \
		os.path.join(PYVISI_TEST_VELOCITY_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_VELOCITY_IMAGES_PATH,\
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestVelocity2DArrowVectorColor(unittest.TestCase, TestVelocity):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_2D))

		self.velocity = Velocity(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, arrow = Arrow.THREE_D, 
				color_mode = ColorMode.VECTOR, lut = Lut.COLOR,
				cell_to_point = False, outline = True)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.velocity

	def testVelocityVectorScale(self):
		self.velocity.setScaleModeByVector()
		self.velocity.setScaleFactor(0.5)
		self.render("TestVelocity2DArrowVectorColor_testVelocityVectorScale.jpg")

class TestVelocity2DArrowScalarColor(unittest.TestCase, TestVelocity):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_2D))

		self.velocity = Velocity(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, arrow = Arrow.TWO_D, 
				color_mode = ColorMode.SCALAR, lut = Lut.COLOR,
				cell_to_point = False, outline = True)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.velocity

	def testScalarScale(self):
		self.velocity.setScaleModeByScalar()
		self.velocity.setScaleFactor(1.5)
		self.render("TestVelocity2DArrowScalarColor_testVelocityScalarScale.jpg")

	def testMask(self):
		self.velocity.setRatio(2)
		self.velocity.randomOn()
		self.render("TestVelocity2DArrowScalarColor_testMask.jpg")

class TestVelocity3DSecondOrder(unittest.TestCase, TestVelocity):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_SECOND_ORDER_3D))

		self.velocity = Velocity(scene = self.scene,
		   		data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, arrow = Arrow.THREE_D, 
				color_mode = ColorMode.VECTOR, lut = Lut.COLOR,
				cell_to_point = False, outline = True)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.velocity

	def testVelocity3DSecondOrder(self):
		self.velocity.setScaleFactor(0.5)
		self.velocity.setRatio(2)
		self.render("TestVelocity3DSecondOrder.jpg")



###############################################################################


if __name__ == '__main__':
    if getMPISizeWorld() == 1:
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity2DArrowVectorColor))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity2DArrowScalarColor))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity3DSecondOrder))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
    else:
        print "run_velocity.py is not executed as more than one processor is used."

