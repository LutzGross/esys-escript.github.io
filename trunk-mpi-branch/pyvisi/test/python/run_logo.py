from esys.pyvisi import Scene, ImageReader, Logo
from esys.pyvisi import LocalPosition
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
PYVISI_TEST_LOGO_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "logo")

MIN_IMAGE_SIZE = 100
LOGO = "access_logo.jpg"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestLogo:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_LOGO_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_LOGO_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestAccessLogo(unittest.TestCase, TestLogo):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
				x_size = X_SIZE, y_size = Y_SIZE)
	
		self.image_reader = ImageReader(ImageFormat.JPG)
		self.image_reader.setImageName(os.path.join(PYVISI_TEST_MESHES_PATH, \
				LOGO))

		self.logo = Logo(scene = self.scene, image_reader = self.image_reader)

	def tearDown(self):
		del self.scene
		del self.image_reader
		del self.logo

	def testImage(self):
		self.logo.setPosition(position = LocalPosition(20,50))
		self.logo.setSize(size = 0.5)
		self.render("TestLogo.jpg")
		

###############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAccessLogo))
	unittest.TextTestRunner(verbosity=2).run(suite)

