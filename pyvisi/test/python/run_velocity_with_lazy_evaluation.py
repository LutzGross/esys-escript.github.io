from esys.pyvisi import DataCollector, Scene, Velocity, VelocityOnPlaneCut
from esys.pyvisi import VelocityOnPlaneClip, Camera
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
PYVISI_TEST_VELOCITY_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "velocity")

MIN_IMAGE_SIZE = 100
FILE_3D = "results.xml"
X_SIZE = 400
Y_SIZE = 400

JPG_RENDERER = Renderer.OFFLINE_JPG
SCALAR_FIELD_POINT_DATA_1 = "scalar1"
SCALAR_FIELD_POINT_DATA_2 = "scalar2"
VECTOR_FIELD_POINT_DATA_1 = "vector"
VECTOR_FIELD_POINT_DATA_2 = "vector2"

class TestVelocityWithLazyEvaluation:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_VELOCITY_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_VELOCITY_IMAGES_PATH,\
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestVelocity(unittest.TestCase, TestVelocityWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testVelocity(self):

		s = Scene(renderer = JPG_RENDERER, num_viewport = 4, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)
		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))
		dc1.setActiveVector(VECTOR_FIELD_POINT_DATA_1)
		dc1.setActiveScalar(SCALAR_FIELD_POINT_DATA_2)

		dc2 = DataCollector(source = Source.XML)
		dc2.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH,\
				FILE_3D))
		dc2.setActiveVector(VECTOR_FIELD_POINT_DATA_1)
		dc2.setActiveScalar(SCALAR_FIELD_POINT_DATA_1)

		dc3 = DataCollector(source = Source.XML)
		dc3.setActiveVector(VECTOR_FIELD_POINT_DATA_2)
		dc3.setActiveScalar(SCALAR_FIELD_POINT_DATA_1)


		# Create a velocity instance in the first viewport.
		v1 = Velocity(scene = s, data_collector = dc1, 
				viewport = Viewport.SOUTH_WEST, 
				arrow = Arrow.THREE_D, color_mode = ColorMode.SCALAR, 
				lut = Lut.COLOR, outline = True)
		v1.setRepresentationToWireframe()
		v1.setScaleFactor(scale_factor = 0.3)
		v1.setScaleModeByScalar()

		# Create a velocity instance in the second viewport.
		v2 = Velocity(scene = s, data_collector = dc2, 
				viewport = Viewport.NORTH_WEST, 
				arrow = Arrow.THREE_D, color_mode = ColorMode.SCALAR, 
				lut = Lut.COLOR, outline = True)
		v2.setScaleModeByScalar()
		v2.setScaleFactor(scale_factor = 0.2)

		# Create a velocity instance in the third viewport.
		v3 = Velocity(scene = s, data_collector = dc2, 
				viewport = Viewport.NORTH_EAST, 
				arrow = Arrow.TWO_D, color_mode = ColorMode.VECTOR, 
				lut = Lut.COLOR, outline = True)

		# Create a velocity instance in the fourth viewport.
		v4 = Velocity(scene = s, data_collector = dc3, 
				viewport = Viewport.SOUTH_EAST,  
				arrow = Arrow.TWO_D, color_mode = ColorMode.SCALAR, 
				lut = Lut.COLOR, outline = True)
		v4.setOpacity(opacity = 0.5)
		v4.setScaleModeByScalar()
		v4.setScaleFactor(scale_factor = 0.2)

		v3.setScaleFactor(scale_factor = 0.2)
		v3.setScaleModeByVector()
		dc3.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))

		self.render("TestVelocityWithLazyEvaluation.jpg")


class TestVelocityOnPlaneCut(unittest.TestCase, TestVelocityWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testVelocityOnPlaneCut(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)
		dc1.setActiveVector(VECTOR_FIELD_POINT_DATA_1)

		# Create a velocity instance.
		vopc1 = VelocityOnPlaneCut(scene = s, data_collector = dc1, 
				arrow = Arrow.THREE_D, color_mode = ColorMode.VECTOR)
		vopc1.setScaleFactor(scale_factor = 0.2)
		vopc1.setPlaneToXY()

		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))

		self.render("TestVelocityOnPlaneCutWithLazyEvaluation.jpg")


class TestVelocityOnPlaneClip(unittest.TestCase, TestVelocityWithLazyEvaluation):
	def tearDown(self):
		del self.scene

	def testVelocityOnPlaneClip(self):
		s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
				y_size = Y_SIZE)
		self.scene = s

		dc1 = DataCollector(source = Source.XML)

		# Create a velocity on plane clip instance.
		vopc1 = VelocityOnPlaneClip(scene = s, data_collector = dc1, 
				arrow = Arrow.THREE_D, color_mode = ColorMode.VECTOR)
		vopc1.setPlaneToXZ()

		c = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c.isometricView()

		vopc1.setScaleFactor(scale_factor = 0.2)
		dc1.setFileName(file_name = os.path.join(PYVISI_TEST_MESHES_PATH, \
				FILE_3D))
		self.render("TestVelocityOnPlaneClipWithLazyEvaluation.jpg")


###############################################################################


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocityOnPlaneCut))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocityOnPlaneClip))
	unittest.TextTestRunner(verbosity=2).run(suite)
