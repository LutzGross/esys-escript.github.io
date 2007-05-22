from esys.pyvisi import DataCollector, Scene, Ellipsoid
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
PYVISI_TEST_ELLIPSOID_REFERENCE_IMAGES_PATH = \
		os.path.join(PYVISI_TEST_DATA_ROOT, "data_reference_images", \
		"ellipsoid")
PYVISI_TEST_ELLIPSOID_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "ellipsoid")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestEllipsoid(unittest.TestCase):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
						x_size = X_SIZE, y_size = Y_SIZE)

		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				os.path.join(PYVISI_TEST_MESHES_PATH, FILE_3D))

		self.ellipsoid = Ellipsoid(scene = self.scene,
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR,
				cell_to_point = False, outline = True)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.ellipsoid

	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_ELLIPSOID_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_ELLIPSOID_IMAGES_PATH,\
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

	def testSphere(self):
		self.ellipsoid.setThetaResolution(5)
		self.ellipsoid.setPhiResolution(10)
		self.render("TestEllipsoid_testSphere.jpg")

	def testTensorGlyph(self):
		self.ellipsoid.setScaleFactor(0.2)
		self.ellipsoid.setMaxScaleFactor(3)
		self.render("TestEllipsoid_testTensorGlyph.jpg")
			

##############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoid))
	unittest.TextTestRunner(verbosity=2).run(suite)

