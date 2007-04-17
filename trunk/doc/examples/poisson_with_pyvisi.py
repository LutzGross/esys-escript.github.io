# $Id: poisson.py 567 2006-02-28 03:58:05Z gross $
from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys.finley import Rectangle
from esys.pyvisi import Scene, DataCollector, Map, Camera
from esys.pyvisi.constant import *

# generate domain:
mydomain = Rectangle(l0=1.,l1=1.,n0=40, n1=20)
# define characteristic function of Gamma^D
x = mydomain.getX()
gammaD = whereZero(x[0])+whereZero(x[1])
# define PDE and get its solution u
mypde = Poisson(domain=mydomain)
mypde.setValue(f=1,q=gammaD)
u = mypde.getSolution()
# write u to an external file
#saveVTK("u.xml",sol=u)

s = Scene(renderer = Renderer.OFFLINE_JPG, x_size = 500, y_size = 500)

dc = DataCollector(source = Source.ESCRIPT)
dc.setData(sol = u)

Map(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST, 
	  lut = Lut.COLOR, cell_to_point = True, outline = True)

c = Camera(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST)

s.render(image_name = "poisson.jpg")
