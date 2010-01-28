########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.pyvisi import DataCollector, Scene, Map, Image, Image, ImageReader
from esys.pyvisi import GlobalPosition
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
PYVISI_TEST_IMAGE_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_sample_images", "image")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
IMAGE = "flinders.jpg"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestImageWithLazyEvaluation:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_IMAGE_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_IMAGE_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestImage(unittest.TestCase, TestImageWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testImage(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create a map instance for the first viewport.
		m1 = Map(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)

		# Create one image reader instance (used in place of data collector).
		ir = ImageReader(ImageFormat.JPG)
			

		# Create one image instance.
		i = Image(scene = s, image_reader = ir)
		i.setOpacity(opacity = 0.9)
		i.translate(0,0,-1.)
		i.setPoint1(GlobalPosition(2,0,0))
		i.setPoint2(GlobalPosition(0,2,0))

		ir.setImageName(image_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				IMAGE)) 
		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))

		self.render("TestImageWithLazyEvaluation.jpg")


###############################################################################


if __name__ == '__main__':
    if getMPISizeWorld() == 1: 
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImage))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
    else:
        print "run_image_with_lazy_evaluation.py is not executed as more than one processor is used."

