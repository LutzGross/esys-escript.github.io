
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
__url__="http://www.uq.edu.au/esscc/escript-finley"

from esys.pyvisi import DataCollector, Scene, Ellipsoid, EllipsoidOnPlaneCut
from esys.pyvisi import EllipsoidOnPlaneClip, Camera
from esys.pyvisi.constant import *
import unittest, os, sys
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
PYVISI_TEST_ELLIPSOID_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "ellipsoid")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestEllipsoidWithLazyEvaluation:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_ELLIPSOID_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_ELLIPSOID_IMAGES_PATH,\
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestEllipsoidScaleResolution(unittest.TestCase, TestEllipsoidWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testEllipsoidScaleResolution(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
        		y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create one ellipsoid instance.
		e1 = Ellipsoid(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST,
				lut = Lut.COLOR, outline = True)
		e1.setScaleFactor(scale_factor = 0.2)
		e1.setPhiResolution(resolution = 20)
		e1.setThetaResolution(resolution = 20)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.isometricView()
		c1.elevation(angle = -20)

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))
		self.render("TestEllipsoidWithLazyEvaluation.jpg")

class TestEllipsoidOnPlaneCut(unittest.TestCase, TestEllipsoidWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testEllipsoidOnPlaneCut(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create one ellipsoid on plance cut instance.
		eopc1 = EllipsoidOnPlaneCut(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
		eopc1.setScaleFactor(scale_factor = 0.1)
		eopc1.setPlaneToXY()
		eopc1.rotateX(angle = 10)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.isometricView()
		c1.elevation(angle = -20)

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))

		self.render("TestEllipsoidOnPlaneCutWithLazyEvaluation.jpg")

class TestEllipsoidOnPlaneClip(unittest.TestCase, TestEllipsoidWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testEllipsoidOnPlaneClip(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s
		
		dc1 = DataCollector(source = Source.XML)

		# Create on ellipsoid on plane clip instance.
		eopc = EllipsoidOnPlaneClip(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
		eopc.setScaleFactor(scale_factor = 0.1)
		eopc.setPlaneToXY()
		eopc.rotateX(angle = 10)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.bottomView()
		c1.azimuth(angle = -90)
		c1.elevation(angle = 20)

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))

		self.render("TestEllipsoidOnPlaneClipWithLazyEvaluation.jpg")


##############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidScaleResolution))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidOnPlaneCut))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidOnPlaneClip))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
