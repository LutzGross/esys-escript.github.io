# Import the necessary modules.
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
from esys.pyvisi import Scene, DataCollector, Map, Camera
from esys.pyvisi.constant import *

PYVISI_EXAMPLE_IMAGES_PATH = "data_sample_images/"
X_SIZE = 400
Y_SIZE = 400
JPG_RENDERER = Renderer.DISPLAY

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
# Create a DataCollector reading directly from escript objects.
dc = DataCollector(source = Source.ESCRIPT)

Map(scene = s, data_collector = dc, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, 
        cell_to_point = False, outline = True)

# ... start iteration:
while t<0.4*100:
      i+=1
      t+=h
      print "time step :",t
      mypde.setValue(Y=qH+rhocp/h*T)
      T=mypde.getSolution()

      dc.setData(temp = T)
      # Create a Map.
      #Map(scene = s, data_collector = dc, 
              #viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, 
              #cell_to_point = False, outline = True)

      # Create a Camera.
      #c= Camera(scene = s, data_collector = dc, 
      #        viewport = Viewport.SOUTH_WEST)
      
      # Render the object.
      s.render(image_name = PYVISI_EXAMPLE_IMAGES_PATH +  
              "diffusion%02d.jpg" % i)

