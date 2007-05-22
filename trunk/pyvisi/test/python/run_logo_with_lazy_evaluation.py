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
PYVISI_TEST_LOGO_REFERENCE_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT, \
		"data_reference_images", "logo")
PYVISI_TEST_LOGO_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "logo")

MIN_IMAGE_SIZE = 100
LOGO = "access_logo.jpg"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestLogoWithLazyEvaluation:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_LOGO_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_LOGO_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestLogo(unittest.TestCase, TestLogoWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testLogo(self):

		# Create a Scene.
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		# Create an ImageReader (in place of DataCollector).
		ir = ImageReader(ImageFormat.JPG)

		# Create a Logo.
		l = Logo(scene = s, image_reader = ir, viewport = Viewport.SOUTH_WEST)
		l.setSize(size = 0.5)

		ir.setImageName(image_name =  os.path.join(PYVISI_TEST_MESHES_PATH, \
				LOGO))
		l.setPosition(position = LocalPosition(50,60))

		# Render the Logo.
		self.render("TestLogoWithLazyEvaluation.jpg")


###########################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLogo))
	unittest.TextTestRunner(verbosity=2).run(suite)

