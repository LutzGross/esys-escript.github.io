
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

from esys.pyvisi import DataCollector, Scene, Contour, ContourOnPlaneCut
from esys.pyvisi import ContourOnPlaneClip, Camera
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
PYVISI_TEST_EXPORTER_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "exporter")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
VRML_RENDERER = Renderer.OFFLINE_VRML
IV_RENDERER = Renderer.OFFLINE_IV

class TestExporter:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_EXPORTER_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join \
				(PYVISI_TEST_EXPORTER_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestVRMLExporter(unittest.TestCase, TestExporter):
	def tearDown(self):
		del self.scene

	def testVRMLExporter(self):
		s = Scene(renderer = VRML_RENDERER, num_viewport = 1, x_size = X_SIZE, 
					y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create one contour instance.
		ctr1 = Contour(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
		ctr1.generateContours(contours = 6, lower_range = 0, upper_range = 0.5)

		cam1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		cam1.elevation(angle = -40)

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))
		self.render("TestVRMLExporter.wrl")

class TestIVExporter(unittest.TestCase, TestExporter):
	def tearDown(self): del self.scene

	def testIVExporter(self):
		s = Scene(renderer = IV_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create one contour on plane cut instance.
		ctropc1 = ContourOnPlaneCut(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
		ctropc1.setPlaneToXY(offset = 0.2)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.elevation(angle = -45)

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))
		ctropc1.generateContours(contours = 8)

		self.render("TestIVExporter.iv")


##############################################################################


if __name__ == '__main__':
    if getMPISizeWorld() == 1: 
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVRMLExporter))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestIVExporter))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
    else:
        print "run_exporter.py is not executed as more than one processor is used."

