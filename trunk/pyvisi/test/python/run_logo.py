from esys.pyvisi import Scene, ImageReader, Logo
from esys.pyvisi import LocalPosition
from esys.pyvisi.constant import *
import unittest, os
from stat import ST_SIZE

PYVISI_TEST_MESHES_PATH = "data_meshes/"
PYVISI_TEST_LOGO_IMAGES_PATH = "data_sample_images/logo/"
MIN_IMAGE_SIZE = 100
LOGO = "access_logo.jpg"

X_SIZE = 400
Y_SIZE = 400

JPG_RENDERER = Renderer.OFFLINE_JPG

class TestLogo:
	def tearDown(self):
		self.scene
		self.image_reader
		self.logo

	def render(self, file):
		self.scene.render(image_name = \
				PYVISI_TEST_LOGO_IMAGES_PATH + file)

		self.failUnless(os.stat(PYVISI_TEST_LOGO_IMAGES_PATH + \
				file)[ST_SIZE] > MIN_IMAGE_SIZE)

class TestAccessLogo(unittest.TestCase, TestLogo):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
				x_size = X_SIZE, y_size = Y_SIZE)
	
		self.image_reader = ImageReader(ImageFormat.JPG)
		self.image_reader.setImageName(PYVISI_TEST_MESHES_PATH + LOGO)

		self.logo = Logo(scene = self.scene, image_reader = self.image_reader)

	def testImage(self):
		self.logo.setPosition(position = LocalPosition(20,50))
		self.logo.setSize(size = 0.5)
		self.render("TestLogo.jpg")
		

###############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAccessLogo))
	unittest.TextTestRunner(verbosity=2).run(suite)

