
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

from esys.pyvisi import DataCollector, Scene, Carpet, Camera
from esys.pyvisi.constant import *
from esys.escript import getMPISizeWorld
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
PYVISI_TEST_CARPET_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_sample_images", "carpet")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestCarpetWithLazyEvaluation:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_CARPET_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_CARPET_IMAGES_PATH,\
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestCarpet(unittest.TestCase, TestCarpetWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testCarpet(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create one carpet instance.
		cpt1 = Carpet(scene = s, data_collector = dc1, 
				warp_mode = WarpMode.SCALAR, lut = Lut.COLOR)
		cpt1.setScaleFactor(0.5)

		c1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c1.isometricView()

		cpt1.setPlaneToXY(0.5)
		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH,\
				FILE_3D))
	
		self.render("TestCarpetWithLazyEvaluation.jpg")


##############################################################################
if __name__ == '__main__':
    if getMPISizeWorld() == 1: 
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCarpet))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)

    else:
        print "run_carpet_with_lazy_evaluation.py is not executed as more than one processor is used."
