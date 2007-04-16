# $Id$
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Brick
from esys.pyvisi import Scene, DataCollector, Map
from esys.pyvisi.constant import *
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
sigma_mises=sqrt(((sigma[0,0]-sigma[1,1])**2+(sigma[1,1]-sigma[2,2])**2+ \
                  (sigma[2,2]-sigma[0,0])**2)/6. \
                   +sigma[0,1]**2 + sigma[1,2]**2 + sigma[2,0]**2)
#... output ...
#saveVTK("deform.xml",disp=u,stress=sigma_mises)
 
s = Scene(renderer = Renderer.ONLINE, x_size = 500, y_size = 500)

dc = DataCollector(source = Source.ESCRIPT)
dc.setData(disp = u,stress = sigma_mises)

Map(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST, 
	  lut = Lut.COLOR, cell_to_point = True, outline = True)

s.render(image_name = "heatedblock.jpg")
