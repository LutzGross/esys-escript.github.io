from esys.pyvisi import Scene
from esys.pyvisi.constant import *
import unittest, os
from stat import ST_SIZE

PYVISI_TEST_SCENE_IMAGES_PATH = "data_sample_images/scene/"
MIN_IMAGE_SIZE = 100
X_SIZE = 400 
Y_SIZE = 400

JPG_RENDERER = Renderer.OFFLINE_JPG

class TestScene:
	def tearDown(self):
		del self.scene

	def render(self, file):
		self.scene.render(image_name = \
				PYVISI_TEST_SCENE_IMAGES_PATH + file)

		self.failUnless(os.stat(PYVISI_TEST_SCENE_IMAGES_PATH + \
				file)[ST_SIZE] > MIN_IMAGE_SIZE)
	
	def setBackground(self, c):
		self.scene.setBackground(color = c)

class TestSceneOneViewport(unittest.TestCase, TestScene):
	def setUp(self): 
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 1, 
				x_size = X_SIZE, y_size = Y_SIZE)
		
	def testRender(self):
		file = "TestSceneOneViewport_testRender.jpg"
		self.render(file)
	
	def testSetBackground(self):
		self.setBackground(Color.GREEN)
		file = "TestSceneOneViewport_testSetBackground.jpg"
		self.render(file)

class TestSceneFourViewports(unittest.TestCase, TestScene):
	def setUp(self):
		self.scene = \
				Scene(renderer = JPG_RENDERER, num_viewport = 4, x_size = 800, 
				y_size = 800) 

	def testRender(self):
		file = "TestSceneFourViewport_testRender.jpg"
		self.render(file)
	

###############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSceneOneViewport))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSceneFourViewports))
	unittest.TextTestRunner(verbosity=2).run(suite)
