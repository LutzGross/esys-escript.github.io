
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

from esys.pyvisi import DataCollector, Scene, Camera, GlobalPosition, Map
from esys.pyvisi.constant import *
from esys.escript import getMPISizeWorld
import unittest, os, gc, sys
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
PYVISI_TEST_CAMERA_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_sample_images", "camera")

MIN_IMAGE_SIZE = 100
FILE_2D = "interior_2D.xml"
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestCamera:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_CAMERA_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_CAMERA_IMAGES_PATH,\
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestCamera2D(unittest.TestCase, TestCamera):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
	    			os.path.join(PYVISI_TEST_MESHES_PATH,  FILE_2D))

	 	self.map = Map(scene = self.scene,
	 			data_collector = self.data_collector,
	 			viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR,
	 			cell_to_point = False, outline = True)

		self.camera = Camera(scene = self.scene, viewport = Viewport.SOUTH_WEST)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.map
		del self.camera

	def test2D(self):
	 	self.camera.azimuth(20)
	 	self.camera.elevation(40)
	 	self.camera.roll(30)
	 	self.camera.dolly(1.5)
		self.render("TestCamera2D_test2D.jpg")

class TestCamera3D(unittest.TestCase, TestCamera):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.map = Map(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR,
				cell_to_point = False, outline = True)

		self.camera = Camera(scene = self.scene,
				viewport = Viewport.SOUTH_WEST)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.map
		del self.camera

	def test3D(self):
		self.camera.setFocalPoint(GlobalPosition(1, 1, 0.3))
		self.camera.setPosition(GlobalPosition(1, 1, 7))
		self.camera.isometricView()
		self.render("TestCamera3D_test3D.jpg")


##############################################################################


if __name__ == '__main__':
        if getMPISizeWorld() == 1: 
	    suite = unittest.TestSuite()
	    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera2D))
	    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera3D))
	    s=unittest.TextTestRunner(verbosity=2).run(suite)
            if not s.wasSuccessful(): sys.exit(1)
        else:
            print "run_camera.py is not executed as more than one processor is used."
