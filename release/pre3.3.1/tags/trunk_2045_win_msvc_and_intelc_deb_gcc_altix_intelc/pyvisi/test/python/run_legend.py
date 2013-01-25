
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

from esys.pyvisi import DataCollector, Scene, Contour, Legend, LocalPosition
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
PYVISI_TEST_LEGEND_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "legend")

MIN_IMAGE_SIZE = 100
FILE_2D = "interior_2D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestLegend(unittest.TestCase):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_2D))

		self.contour = Contour(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR,
				cell_to_point = False, outline = True)
			

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.contour
		del self.legend

	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_LEGEND_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_LEGEND_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

	def testScalarLegend(self):
		self.legend = Legend(scene = self.scene, \
				data_collector = self.data_collector, \
				viewport = Viewport.SOUTH_WEST,\
		        lut = Lut.COLOR, legend = LegendType.SCALAR)

		self.legend.setOrientationToHorizontal()
		self.legend.setTitle(title = "Scalar Bar")
		self.legend.setPosition(LocalPosition(50, 5))

		self.render("TestLegend_testScalarLegend.jpg")

	def testVectorLegend(self):
		self.legend = Legend(scene = self.scene, \
				data_collector = self.data_collector, \
				viewport = Viewport.SOUTH_WEST,\
		        lut = Lut.COLOR, legend = LegendType.VECTOR)

		self.legend.setOrientationToVertical()
		self.legend.setTitle(title = "Scalar Bar")
		self.legend.setPosition(LocalPosition(5, 50))

		self.render("TestLegend_testVectorLegend.jpg")



##############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLegend))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)

