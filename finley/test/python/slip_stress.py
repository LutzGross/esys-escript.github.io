# $Id$
"""
calculation of the stress distribution around a fault from the slip on the fault

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, Louise Kettle"
__copyright__="""  Copyright (c) 2006, 2007 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"

from esys.escript import *
from esys.escript.pdetools import SaddlePointProblem
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Brick
from esys.pyvisi import Scene, DataCollector, Contour, Camera, Velocity, Text2D, LocalPosition
from esys.pyvisi.constant import *
import os

JPG_RENDERER = Renderer.OFFLINE_JPG
ne = 20          # number of elements in spatial direction

def lockToGrid(x,l,n):
   h=l/n
   return int(x/h+0.5)*h


width  = 100000.
height =  30000.

ne_w=int((ne/height)*width+0.5)
rho=0.
lmbd=1.7e11
mu=1.7e11
g=9.81

# fault location:

fstart =  [lockToGrid(50000.0,width,ne_w), lockToGrid(40000.0,width,ne_w), lockToGrid(8000.,height,ne)]
fend =  [lockToGrid(50000.0,width,ne_w), lockToGrid(60000.0,width,ne_w), lockToGrid(20000.,height,ne)]
s=[0.,1.,0.]

print "=== generate mesh over %s x %s x %s ==="%(width,width,height)
dom=Brick(l0=width, l1=width, l2=height, n0=ne_w, n1=ne_w, n2=ne)
print "  total number of elements = ",ne*ne_w*ne_w

print "=== prepare PDE coefficients ==="
# fixed displacements:
d=dom.getDim()
x=dom.getX()
mask=whereZero(x[d-1]-inf(x[d-1]))*numarray.ones((d,))

# map the fault into the unit square:
x_hat=x-fstart
p=[0.,0.,0.]
q=1.
for i in range(d):
     p[i]=fend[i]-fstart[i]
     if abs(p[i])>0: 
         x_hat[i]/=p[i]
         q=whereNonNegative(x_hat[i])*whereNonPositive(x_hat[i]-1.)*x_hat[i]*(1.-x_hat[i])*4.*q
     else:
         flip=i
         q=whereZero(x_hat[i],1.e-3*Lsup(x_hat[i]))*q
g_s=grad(q*s)
sigma0=(mu*symmetric(g_s)+lmbd*trace(g_s)*kronecker(d))*(whereNegative(interpolate(x_hat[flip],g_s.getFunctionSpace()))-0.5)
pde=LinearPDE(dom)
pde.setSymmetryOn()
A =Tensor4(0.,Function(dom))
for i in range(d):
   for j in range(d):
      A[i,j,j,i] += mu
      A[i,j,i,j] += mu
      A[i,i,j,j] += lmbd
pde.setValue(A=A,q=mask,Y=-kronecker(Function(dom))[d-1]*g*rho,X=sigma0)

print "=== solve pde === "
u=pde.getSolution(verbose=True)

print "=== calculate stress ==="
g_s=grad(u)
sigma=mu*symmetric(g_s)+lmbd*trace(g_s)*kronecker(d)

print "=== start rendering ==="
# Create a Scene.
s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = 800, y_size = 800)
dc1 = DataCollector(source = Source.ESCRIPT)
dc1.setData(disp=u,sigma=sigma,cfs=sigma[0,1]-0.4*sigma[0,0])

# Create a Contour.
ctr1 = Contour(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST,
        lut = Lut.COLOR, cell_to_point = True, outline = True)
ctr1.generateContours(contours = 10)
ctr1.setOpacity(0.5)

# Create velocity
vopc1 = Velocity(scene = s, data_collector = dc1,
        color_mode = ColorMode.VECTOR,
        arrow = Arrow.THREE_D, lut = Lut.COLOR, cell_to_point = False,
        outline = False) 
vopc1.setScaleFactor(scale_factor = 10000.)
vopc1.setRatio(2)
vopc1.randomOn()

# Create a Camera.
cam1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
cam1.elevation(angle = -50)

# add some text
t2 = Text2D(scene = s, text = "CFS and displacement around vertical fault")
t2.setPosition(LocalPosition(200,30))
t2.setColor(color = Color.BLACK)
t2.setFontSize(size = 20)
t2.setFontToArial()
t2.shadowOn()

# Render the object.
s.render("stress_contour.jpg")
print "=== finished rendering ==="
