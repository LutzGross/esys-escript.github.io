# $Id: diffusion.py 575 2006-03-03 03:33:07Z lkettle $
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
from esys.pyvisi import Scene, DataCollector, Map, Camera
from esys.pyvisi.constant import *

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

s = Scene(renderer = Renderer.OFFLINE_JPG, x_size = 500, y_size = 500)
dc = DataCollector(source = Source.ESCRIPT)

# ... start iteration:
while t<tend:
      i+=1
      t+=h
      print "time step :",t
      mypde.setValue(Y=qH+rhocp/h*T)
      T=mypde.getSolution()
      #saveVTK("T.%d.xml"%i,temp=T)

      dc.setData(temp = T)
      Map(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST, 
              lut = Lut.COLOR, cell_to_point = False, outline = True)

      c = Camera(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST)

      s.render(image_name = "diffusion_%02d.jpg" % i)

