# Import the necessary modules.
from esys.escript import *
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle, Brick
from numarray import identity,zeros,ones
from esys.pyvisi import Scene, DataCollector, Map, Velocity, Ellipsoid, Camera
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

PYVISI_TEST_ESCRIPT_REFERENCE_IMAGES_PATH = os.path.join(PYVISI_TEST_DATA_ROOT,\
		"data_reference_images", "escript")
PYVISI_TEST_ESCRIPT_IMAGES_PATH = os.path.join(PYVISI_WORKDIR, \
		"data_sample_images", "escript")

MIN_IMAGE_SIZE = 100
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.OFFLINE_JPG

class TestEscript:
	def render(self, file):
		self.scene.render(image_name = \
				os.path.join(PYVISI_TEST_ESCRIPT_IMAGES_PATH, file))

		self.failUnless(os.stat(os.path.join(PYVISI_TEST_ESCRIPT_IMAGES_PATH, \
				file))[ST_SIZE] > MIN_IMAGE_SIZE)

class TestEscriptMap(unittest.TestCase, TestEscript):
	def tearDown(self):
		del self.scene

	def testEscriptMap(self):
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
		s = Scene(renderer = JPG_RENDERER, x_size = X_SIZE, y_size = Y_SIZE)
		self.scene = s

		# Create a DataCollector reading directly from escript objects.
		dc = DataCollector(source = Source.ESCRIPT)

		# Create a Map.
		m = Map(scene = s, data_collector = dc, \
				viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, \
				cell_to_point = False, outline = True)

		# Create a Camera.
		c = Camera(scene = s, viewport = Viewport.SOUTH_WEST)

		# ... start iteration:
		while t<0.4:
			i+=1
			t+=h
			mypde.setValue(Y=qH+rhocp/h*T)
			T=mypde.getSolution()

			dc.setData(temp = T)

			# Render the object.
			self.render("TestEscriptMap%02d.jpg" % i)

class TestEscriptVelocity(unittest.TestCase, TestEscript):
	def tearDown(self):
		del self.scene

	def testEscriptVelocity(self):
		#... set some parameters ...
		lam=1.
		mu=0.1
		alpha=1.e-6
		xc=[0.3,0.3,1.]
		beta=8.
		T_ref=0.
		T_0=1.

		#... generate domain ...
		mydomain = Brick(l0=1.,l1=1., l2=1.,n0=10, n1=10, n2=10)
		x=mydomain.getX()
		#... set temperature ...
		T=T_0*exp(-beta*length(x-xc))
		#... open symmetric PDE ...
		mypde=LinearPDE(mydomain)
		mypde.setSymmetryOn()

		#... set coefficients ...
		C=Tensor4(0.,Function(mydomain))
		for i in range(mydomain.getDim()):
			for j in range(mydomain.getDim()):
				C[i,i,j,j]+=lam
				C[j,i,j,i]+=mu
				C[j,i,i,j]+=mu

		msk=whereZero(x[0])*[1.,0.,0.] \
				+whereZero(x[1])*[0.,1.,0.] \
				+whereZero(x[2])*[0.,0.,1.]
		sigma0=(lam+2./3.*mu)*alpha*(T-T_ref)*kronecker(mydomain)
		mypde.setValue(A=C,X=sigma0,q=msk)

		#... solve pde ...
		u=mypde.getSolution()
		#... calculate von-Misses
		g=grad(u)
		sigma=mu*(g+transpose(g))+lam*trace(g)*kronecker(mydomain)-sigma0
		sigma_mises=sqrt(((sigma[0,0]-sigma[1,1])**2+\
				(sigma[1,1]-sigma[2,2])**2+ \
				(sigma[2,2]-sigma[0,0])**2)/6. \
				+sigma[0,1]**2 + sigma[1,2]**2 + sigma[2,0]**2)

		# Create a Scene.
		s = Scene(renderer = JPG_RENDERER, x_size = X_SIZE, y_size = Y_SIZE)
		self.scene = s

		# Create a DataCollector reading directly from an escript object.
		dc = DataCollector(source = Source.ESCRIPT)
		dc.setData(disp = u, stress = sigma_mises)

		# Create a Velocity.
		v = Velocity(scene = s, data_collector = dc, \
				viewport = Viewport.SOUTH_WEST, \
				arrow = Arrow.THREE_D, color_mode = ColorMode.SCALAR, \
				lut = Lut.COLOR, cell_to_point = True, outline = True)
		v.setScaleFactor(scale_factor = 0.3)

		# Create a Camera.
		c = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c.isometricView()

		# Render the object.
		self.render("TestEscriptVelocity.jpg")

ne=32 # number of cells in x_0 and x_1 directions
width=10000. # length in x_0 and x_1 directions
lam=3.462e9
mu=3.462e9
rho=1154.
tend=60
h=(1./5.)*sqrt(rho/(lam+2*mu))*(width/ne)

U0=0.01 # amplitude of point source
class TestEscriptEllipsoid(unittest.TestCase, TestEscript):
	def tearDown(self):
		del self.scene

	def testEscriptEllipsoid(self):
		mydomain=Brick(ne,ne,10,l0=width,l1=width,l2=10.*width/32.)
		self.wavePropagation(mydomain,h,tend,lam,mu,rho,U0)

	def wavePropagation(self,domain,h,tend,lam,mu,rho,U0):
		x=domain.getX()
		# ... open new PDE ...
		mypde=LinearPDE(domain)
		mypde.setSolverMethod(LinearPDE.LUMPING)
		kronecker=identity(mypde.getDim())

		#  spherical source at middle of bottom face

		xc=[width/2.,width/2.,0.]
		# define small radius around point xc
		# Lsup(x) returns the maximum value of the argument x
		src_radius = 0.1*Lsup(domain.getSize())
		dunit=numarray.array([1.,0.,0.]) # defines direction of point source

		mypde.setValue(D=kronecker*rho)
		# ... set initial values ....
		n=0
		# initial value of displacement at point source is constant (U0=0.01)
		# for first two time steps
		u=U0*whereNegative(length(x-xc)-src_radius)*dunit
		u_last=U0*whereNegative(length(x-xc)-src_radius)*dunit
		t=0

		# define the location of the point source 
		L=Locator(domain,numarray.array(xc))
		# find potential at point source
		u_pc=L.getValue(u)

		u_pc_x = u_pc[0]
		u_pc_y = u_pc[1]
		u_pc_z = u_pc[2]

		# open file to save displacement at point source
		#u_pc_data=open('./data/U_pc.out','w')
		#u_pc_data.write("%f %f %f %f\n"%(t,u_pc_x,u_pc_y,u_pc_z))

		# Create a Scene.
		s = Scene(renderer = JPG_RENDERER, x_size = X_SIZE, y_size = Y_SIZE)
		self.scene = s

		# Create a DataCollector reading directly from escript objects.
		dc = DataCollector(source = Source.ESCRIPT)

		# Create an Ellipsoid.
		e = Ellipsoid(scene = s, data_collector = dc, 
		viewport = Viewport.SOUTH_WEST,
		lut = Lut.COLOR, cell_to_point = True, outline = True)
		e.setScaleFactor(scale_factor = 0.7)
		e.setMaxScaleFactor(max_scale_factor = 1000)
		e.setRatio(ratio = 10)

		# Create a Camera.
		c = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
		c.isometricView()

		while t<0.5:
			# ... get current stress ....
			g=grad(u)
			stress=lam*trace(g)*kronecker+mu*(g+transpose(g))
			# ... get new acceleration ....
			mypde.setValue(X=-stress)          
			a=mypde.getSolution()
			# ... get new displacement ...
			u_new=2*u-u_last+h**2*a
			# ... shift displacements ....
			u_last=u
			u=u_new
			t+=h
			n+=1
			u_pc=L.getValue(u)

			u_pc_x=u_pc[0]
			u_pc_y=u_pc[1]
			u_pc_z=u_pc[2]

			# ... save current acceleration in units of gravity and 
			# displacements 
			if n==1 or n%10==0: 
				dc.setData(acceleration = length(a)/9.81, displacement = u, 
						tensor = stress, Ux = u[0])

				# Render the object.
				self.render("TestEscriptEllipsoid%02d.jpg" % (n/10))


if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptMap))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptVelocity))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptEllipsoid))
	unittest.TextTestRunner(verbosity=2).run(suite)

