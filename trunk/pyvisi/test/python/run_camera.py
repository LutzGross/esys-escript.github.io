from esys.pyvisi import DataCollector, Scene, Camera, GlobalPosition, Map
from esys.pyvisi.constant import *
import unittest, os
from stat import ST_SIZE

<<<<<<< .mine
try:
	PYVISI_WORKDIR=os.environ['PYVISI_WORKDIR']
except KeyError:
	PYVISI_WORKDIR='.'
try:
	PYVISI_TEST_DATA_ROOT=os.environ['PYVISI_TEST_DATA_ROOT']
except KeyError:
	PYVISI_TEST_DATA_ROOT='.'

PYVISI_TEST_MESHES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, "data_meshes")
PYVISI_TEST_CAMERA_REFERENCE_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_reference_images", "camera")
PYVISI_TEST_CAMERA_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "camera")

=======
try:
     PYVISI_WORKDIR=os.environ['PYVISI_WORKDIR']
except KeyError:
     PYVISI_WORKDIR='.'
     
PYVISI_TEST_MESHES_PATH = s.path.join(PYVISI_WORKDIR,"data_meshes")
PYVISI_TEST_CAMERA_IMAGES_PATH = os.path.join(PYVISI_WORKDIR,"data_sample_images","camera")
>>>>>>> .r1142
MIN_IMAGE_SIZE = 100
FILE_2D = "interior_2D.xml"
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestCamera:
	def tearDown(self):
		self.scene
		self.data_collector
		self.map
		self.camera

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

	def test3D(self):
		self.camera.setFocalPoint(GlobalPosition(1, 1, 0.3))
		self.camera.setPosition(GlobalPosition(1, 1, 7))
		self.camera.isometricView()
		self.render("TestCamera3D_test3D.jpg")


##############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera2D))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera3D))
	unittest.TextTestRunner(verbosity=2).run(suite)

