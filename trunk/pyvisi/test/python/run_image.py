
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
PYVISI_TEST_IMAGE_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "image")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
IMAGE = "flinders.jpg"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestImage:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_IMAGE_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_IMAGE_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestImageOnAMap(unittest.TestCase, TestImage):
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

		self.image_reader = ImageReader(ImageFormat.JPG)
		self.image_reader.setImageName(os.path.join(PYVISI_TEST_MESHES_PATH,\
				IMAGE))

		self.image = Image(scene = self.scene, image_reader = self.image_reader)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.map
		del self.image
		del self.image_reader

	def testImage(self):
		self.image.translate(0, 0, -1)
		self.image.setPoint1(GlobalPosition(2, 0, 0))
		self.image.setPoint2(GlobalPosition(0, 2, 0))
		self.image.setOpacity(0.6)
		self.render("TestImage.jpg")
		

###############################################################################


if __name__ == '__main__':
    if getMPISizeWorld() == 1: 
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImageOnAMap))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
    else:
        print "run_image.py is not executed as more than one processor is used."

