from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle, Brick
from esys.pyvisi import Scene, DataCollector, Map, Camera, Velocity, Ellipsoid
from esys.pyvisi.constant import *
import unittest, os
from stat import ST_SIZE

PYVISI_TEST_ESCRIPT_IMAGES_PATH = "data_sample_images/escript/"
MIN_IMAGE_SIZE = 100

X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestReadingEscriptObject:
	def tearDown(self):
		self.s

	def render(self, file):
		self.s.render(image_name = PYVISI_TEST_ESCRIPT_IMAGES_PATH + file)
			
		self.failUnless(os.stat(PYVISI_TEST_ESCRIPT_IMAGES_PATH + \
				file)[ST_SIZE] > MIN_IMAGE_SIZE)

class TestEscriptWithPointData(unittest.TestCase, TestReadingEscriptObject):
	def testPointData(self):
		#... set some parameters ...
		xc=[0.02,0.002]
		r=0.001
		qc=50.e6
		Tref=0.
		rhocp=2.6e6
		eta=75.
		kappa=240.
		tend=5.
		# ... time, time step size and counter ...
		t=0
		h=0.1
		i=0
		#... generate domain ...
		mydomain = Rectangle(l0=0.05,l1=0.01,n0=250, n1=50)
		#... open PDE ...
		mypde=LinearPDE(mydomain)
		mypde.setSymmetryOn()
		mypde.setValue(A=kappa*kronecker(mydomain),D=rhocp/h,d=eta,y=eta*Tref)
		# ... set heat source: ....
		x=mydomain.getX()
		qH=qc*whereNegative(length(x-xc)-r)
		# ... set initial temperature ....
		T=Tref

		# Create a Scene.
		self.s  = Scene(renderer = JPG_RENDERER, x_size = X_SIZE, \
				y_size = Y_SIZE)
		s  = self.s
		# Create a DataCollector reading directly from escript objects.
		dc = DataCollector(source = Source.ESCRIPT)

		# ... start iteration:
		while t<0.4:
			  i+=1
			  t+=h
			  print "time step :",t
			  mypde.setValue(Y=qH+rhocp/h*T)
			  T=mypde.getSolution()

			  dc.setData(temp = T)
			  # Create a Map.
			  Map(scene = s, data_collector = dc, 
					  viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, 
					  cell_to_point = False, outline = True)

			  # Create a Camera.
			  c= Camera(scene = s, data_collector = dc, 
					  viewport = Viewport.SOUTH_WEST)
			  
			  # Render the object.
			  self.render("diffusion_%02d.jpg" % i)


#class TestEscriptWithPointAndCellData(unittest.TestCase, TestReadingEscriptObject):
#	def testPointAndCellData(self):
		
###############################################################################
if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptWithPointData))
	unittest.TextTestRunner(verbosity=2).run(suite)
