
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

from esys.pyvisi import Scene, DataCollector, Map, Camera, Velocity, Legend 
from esys.pyvisi import Movie, LocalPosition
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
PYVISI_TEST_MOVIE_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "movie")

MIN_IMAGE_SIZE = 100
X_SIZE = 400
Y_SIZE = 400

SCALAR_FIELD_POINT_DATA = "temp"
FILE_2D = "tempvel-"
IMAGE_NAME = "TestMovie"
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestMovie:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_MOVIE_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_MOVIE_IMAGES_PATH,\
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestGenerateMovie(unittest.TestCase, TestMovie):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setActiveScalar(scalar = SCALAR_FIELD_POINT_DATA)


		# Create a Contour.
		self.map = Map(scene = self.scene, 
				data_collector = self.data_collector, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, 
				cell_to_point = False, outline = True)

		self.velocity = Velocity(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, arrow = Arrow.TWO_D, 
				color_mode = ColorMode.VECTOR, lut = Lut.COLOR,
				cell_to_point = False, outline = True)

		self.velocity.setScaleFactor(scale_factor = 0.07)
		self.velocity.setRatio(ratio = 8)
		self.velocity.setColor(color = Color.BLACK)

		self.mov = Movie()
		self.lst = []
		self.cam = Camera(scene = self.scene, viewport = Viewport.SOUTH_WEST)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.velocity
		del self.map
		del self.cam
		del self.mov
		del self.lst

	def testMovieRange(self):

		for i in range(938, 949):
			self.data_collector.setFileName(file_name = \
					os.path.join(PYVISI_TEST_MESHES_PATH, FILE_2D + "%06d.vtu")
					% i)
			
			self.render(IMAGE_NAME + "%06d.jpg" % i)

		self.mov.imageRange(input_directory = PYVISI_TEST_MOVIE_IMAGES_PATH, 
				first_image = IMAGE_NAME + "000938.jpg", 
				last_image = IMAGE_NAME + "000948.jpg")

		self.mov.makeMovie(os.path.join(PYVISI_TEST_MOVIE_IMAGES_PATH, \
				"movie_testMovieRange.mpg"))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_MOVIE_IMAGES_PATH,\
				"movie_testMovieRange.mpg"))[ST_SIZE] > MIN_IMAGE_SIZE)


	def testMovieList(self):

		for i in range(938, 949):
			self.data_collector.setFileName(file_name = \
					os.path.join(PYVISI_TEST_MESHES_PATH, FILE_2D + "%06d.vtu")
					% i)
			
			self.render(IMAGE_NAME + "%06d.jpg" % i)
			self.lst.append(IMAGE_NAME + "%06d.jpg" % i)

		self.mov.imageList(input_directory = PYVISI_TEST_MOVIE_IMAGES_PATH, 
				image_list = self.lst)

		self.mov.makeMovie(os.path.join(PYVISI_TEST_MOVIE_IMAGES_PATH, \
				"movie_testMovieList.mpg"))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_MOVIE_IMAGES_PATH,\
				"movie_testMovieList.mpg"))[ST_SIZE] > MIN_IMAGE_SIZE)


###############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGenerateMovie))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)

