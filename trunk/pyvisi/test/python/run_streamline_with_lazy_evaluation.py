from esys.pyvisi import DataCollector, Scene, StreamLine, GlobalPosition
from esys.pyvisi import Camera
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
PYVISI_TEST_STREAMLINE_REFERENCE_IMAGES_PATH = \
		os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_reference_images", "streamline")
PYVISI_TEST_STREAMLINE_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "streamline")

MIN_IMAGE_SIZE = 100
FILE_3D = "results.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestStreamLineWithLazyEvaluation:
	def tearDown(self):
		self.scene

	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_STREAMLINE_IMAGES_PATH, file))

		self.failUnless(os.stat(\
				os.path.join(PYVISI_TEST_STREAMLINE_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestStreamLine(unittest.TestCase, TestStreamLineWithLazyEvaluation):
	def testStreamLine(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create one streamline instance for the first viewport.
		sl1 = StreamLine(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, 
				outline = True, color_mode = ColorMode.VECTOR)
		sl1.setTubeRadius(radius = 0.01)

		cam1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		cam1.elevation(angle = -40)
		dc1.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.render("TestStreamLineWithLazyEvaluation.jpg")


##############################################################################
if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLine))
	unittest.TextTestRunner(verbosity=2).run(suite)

