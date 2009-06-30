
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Lutz Gross, l.gross@uq.edu.au
Author: John Ngui, john.ngui@uq.edu.au
"""

# Import the necessary modules.
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Brick
from numpy import identity,zeros,ones
from esys.pyvisi import Scene, DataCollector, Ellipsoid, Camera
from esys.pyvisi.constant import *
import os

PYVISI_EXAMPLE_IMAGES_PATH = "images_out"
if not os.path.isdir(PYVISI_EXAMPLE_IMAGES_PATH) and getMPIRankWorld()==0: os.mkdir(PYVISI_EXAMPLE_IMAGES_PATH)

X_SIZE = 400
Y_SIZE = 300
JPG_RENDERER = Renderer.OFFLINE_JPG # change to Renderer.ONLINE_JPG to interact with visualiztion window

ne=32          # number of cells in x_0 and x_1 directions
width=10000.  # length in x_0 and x_1 directions
lam=3.462e9
mu=3.462e9
rho=1154.
tend=0.5     # to ran a full simulation change tend to 60.
alpha=0.3

U0=0.01 # amplitude of point source

def wavePropagation(domain,h,tend,lam,mu,rho,U0):
   x=domain.getX()
   # ... open new PDE ...
   mypde=LinearPDE(domain)
   mypde.getSolverOptions().setSolverMethod(mypde.getSolverOptions().LUMPING)
   kronecker=identity(mypde.getDim())

   #  spherical source at middle of bottom face
   xc=[width/2.,width/2.,0.]
   # define small radius around point xc
   src_radius = 0.03*width
   print "src_radius = ",src_radius

   dunit=numpy.array([1.,0.,0.]) # defines direction of point source

   mypde.setValue(D=kronecker*rho, q=whereNegative(length(x-xc)-src_radius)*dunit)
   # ... set initial values ....
   n=0
   # initial value of displacement at point source is constant (U0=0.01)
   # for first two time steps
   u=Vector(0.,Solution(domain))
   u_last=Vector(0.,Solution(domain))
   t=0

   # Create a Scene.
   s = Scene(renderer = JPG_RENDERER, x_size = X_SIZE, y_size = Y_SIZE)

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

   while t<0.4:
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
   while t<tend:
     t+=h
     # ... get current stress ....
     g=grad(u)
     stress=lam*trace(g)*kronecker+mu*(g+transpose(g))
     # ... get new acceleration ....
     amplitude=U0*2*exp(1)/alpha**2*(1-5*(t/alpha)**2+2*(t/alpha)**4)*exp(-(t/alpha)**2)
     mypde.setValue(X=-stress, r=dunit*amplitude)
     a=mypde.getSolution()
     # ... get new displacement ...
     u_new=2*u-u_last+h**2*a
     # ... shift displacements ....
     u_last=u
     u=u_new
     n+=1
     print n,"-th time step t ",t
     # ... save current acceleration in units of gravity and displacements 
     if n==1 or n%10==0: 
         dc.setData(acceleration = length(a)/9.81, displacement = u, 
                 tensor = stress, Ux = u[0])
        
         # Render the object.
         s.render(image_name = os.path.join(PYVISI_EXAMPLE_IMAGES_PATH, "wave%02d.jpg") % (n/10))

mydomain=Brick(ne,ne,10,l0=width,l1=width,l2=10.*width/ne)
h=inf(1./5.)*inf(sqrt(rho/(lam+2*mu))*mydomain.getSize())
print "time step size = ",h
wavePropagation(mydomain,h,tend,lam,mu,rho,U0)

