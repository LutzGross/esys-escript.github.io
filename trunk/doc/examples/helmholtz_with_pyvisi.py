# $Id: helmholtz.py 575 2006-03-03 03:33:07Z lkettle $
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
from esys.pyvisi import Scene, DataCollector, Map, Camera
from esys.pyvisi.constant import *

#... set some parameters ...
kappa=1.
omega=0.1
eta=10.
#... generate domain ...
mydomain = Rectangle(l0=5.,l1=1.,n0=50, n1=10)
#... open PDE and set coefficients ...
mypde=LinearPDE(mydomain)
mypde.setSymmetryOn()
n=mydomain.getNormal()
x=mydomain.getX()
mypde.setValue(A=kappa*kronecker(mydomain),D=omega,Y=omega*x[0], \
               d=eta,y=kappa*n[0]+eta*x[0])
#... calculate error of the PDE solution ...
u=mypde.getSolution()
print "error is ",Lsup(u-x[0])
# output should be similar to "error is 1.e-7"
#saveVTK("x0.xml",sol=u)
 
s = Scene(renderer = Renderer.OFFLINE_JPG, x_size = 800, y_size = 600)

dc = DataCollector(source = Source.ESCRIPT)
dc.setData(sol = u)

Map(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST, 
	  lut = Lut.COLOR, cell_to_point = False, outline = True)

c = Camera(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST)

s.render(image_name = "helmholtz.jpg")
