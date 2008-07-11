from esys.pyvisi import DataCollector, Scene, Map, Rectangle
from esys.pyvisi import GlobalPosition
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
PYVISI_TEST_IMAGE_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "rectangle")

MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestRectangle:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_IMAGE_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_IMAGE_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestRectangleOnAMap(unittest.TestCase, TestRectangle):
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
		self.map.setOpacity(0.2)

		self.rectangle = Rectangle(scene = self.scene, viewport = Viewport.SOUTH_WEST)

	def tearDown(self):
		del self.scene
		del self.data_collector
		del self.map
		del self.rectangle

	def testRectangle(self):
		self.rectangle.setCenter(GlobalPosition(1.5,2,0.8))
		self.rectangle.setXLength(2.5)
		self.rectangle.setYLength(1)
		self.rectangle.setZLength(0.3)
		self.rectangle.setColor(Color.GREY)

		self.render("TestRectangle.jpg")
		

###############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRectangleOnAMap))
	s=unittest.TextTestRunner(verbosity=2).run(suite)
        if not s.wasSuccessful(): sys.exit(1)
