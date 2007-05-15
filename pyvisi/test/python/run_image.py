from esys.pyvisi import DataCollector, Scene, Map, Image, Image, ImageReader
from esys.pyvisi import GlobalPosition
from esys.pyvisi.constant import *
import unittest, os
from stat import ST_SIZE

PYVISI_TEST_MESHES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT,"data_meshes")
PYVISI_TEST_IMAGE_IMAGES_PATH = "data_sample_images/image/"
MIN_IMAGE_SIZE = 100
FILE_3D = "interior_3D.xml"
IMAGE = "flinders.jpg"

X_SIZE = 400
Y_SIZE = 400

JPG_RENDERER = Renderer.OFFLINE_JPG

class TestImage:
	def tearDown(self):
		self.scene
		self.data_collector
		self.map
		self.image
		self.image_reader

	def render(self, file):
		self.scene.render(image_name = \
				PYVISI_TEST_IMAGE_IMAGES_PATH + file)

		self.failUnless(os.stat(PYVISI_TEST_IMAGE_IMAGES_PATH + \
				file)[ST_SIZE] > MIN_IMAGE_SIZE)

class TestImageOnAMap(unittest.TestCase, TestImage):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1,
				x_size = X_SIZE, y_size = Y_SIZE)
	
		self.data_collector = DataCollector(source = Source.XML)
		self.data_collector.setFileName(file_name = \
				PYVISI_TEST_MESHES_PATH + FILE_3D)
		
		self.map = Map(scene = self.scene, 
				data_collector = self.data_collector,
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR,
				cell_to_point = False, outline = True)

		self.image_reader = ImageReader(ImageFormat.JPG)
		self.image_reader.setImageName(PYVISI_TEST_MESHES_PATH + IMAGE)

		self.image = Image(scene = self.scene, image_reader = self.image_reader)

	def testImage(self):
		self.image.translate(0, 0, -1)
		self.image.setPoint1(GlobalPosition(2, 0, 0))
		self.image.setPoint2(GlobalPosition(0, 2, 0))
		self.image.setOpacity(0.6)
		self.render("TestImage.jpg")
		

###############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImageOnAMap))
	unittest.TextTestRunner(verbosity=2).run(suite)

