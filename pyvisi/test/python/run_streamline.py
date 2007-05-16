from esys.pyvisi import DataCollector, Scene, StreamLine, GlobalPosition
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
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestStreamLine:
	def tearDown(self):
		self.scene
		self.data_collector
		self.streamline

	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_STREAMLINE_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(\
				PYVISI_TEST_STREAMLINE_IMAGES_PATH, file))[ST_SIZE] > \
				MIN_IMAGE_SIZE)

class TestStreamLinePointSource(unittest.TestCase, TestStreamLine):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.streamline = StreamLine(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.VECTOR,
				lut = Lut.COLOR, cell_to_point = False, outline = True)

	def testPointSource(self):
		self.streamline.setPointSourceRadius(0.3)
		self.streamline.setPointSourceCenter(GlobalPosition(1.3, 1.3, 0.4))
		self.streamline.setPointSourceNumberOfPoints(5)
		self.render("TestStreamLinePointSource.jpg")

class TestStreamLineModule(unittest.TestCase, TestStreamLine):	
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.streamline = StreamLine(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.VECTOR,
				lut = Lut.COLOR, cell_to_point = False, outline = True)

	def testStreamLineModule(self):
		self.streamline.setMaximumPropagationTime(20)
		self.streamline.setStepLength(0.1)
		self.streamline.setIntegrationStepLength(0.1)
		self.streamline.setIntegrationToBothDirections()
		self.render("TestStreamLineModule.jpg")

class TestStreamLineTube(unittest.TestCase, TestStreamLine):	
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.streamline = StreamLine(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.VECTOR,
				lut = Lut.COLOR, cell_to_point = False, outline = True)

	def testSetTubeRadiusToVaryByVector(self):
		self.streamline.setTubeRadius(0.05)
		self.streamline.setTubeNumberOfSides(3)
		self.streamline.setTubeRadiusToVaryByVector()
		self.render("TestStreamLineTube_testSetTubeRadiusToVaryByVector.jpg")

	def testSetTubeRadiusToVaryBeScalar(self):
		self.streamline.setTubeRadiusToVaryByScalar()
		self.render("TestStreamLineTube_testSetTubeRadiusToVaryByScalar.jpg")
		

##############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLinePointSource))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineModule))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineTube))
	unittest.TextTestRunner(verbosity=2).run(suite)

