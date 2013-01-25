# $Id$

"""
calculation of the stress distribution around a fault from the slip on the fault

e.g. use slip_stress_mesh.py to generate mesh

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, Louise Kettle"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
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

def lockToGrid(x,l,n):
   h=l/n
   print x,l,n,"->",h
   return int(x/h+0.5)*h


ne = 15
width  = 100000.
height =  30000.

ne_w=int((ne/height)*width+0.5)
rho=0.
lmbd=1.7e11
mu=1.7e11
g=9.81
print "total ne = ",ne*ne_w*ne_w

# fault:

fstart =  [lockToGrid(50000.0,width,ne_w), lockToGrid(40000.0,width,ne_w), lockToGrid(8000.,height,ne)]
fend =  [lockToGrid(50000.0,width,ne_w), lockToGrid(60000.0,width,ne_w), lockToGrid(20000.,height,ne)]
s=[0.,1.,0.]

# fstart =  [lockToGrid(30000.0,width,ne_w), lockToGrid(30000.0,width,ne_w), lockToGrid(20000.,height,ne)]
# fend =  [lockToGrid(70000.0,width,ne_w), lockToGrid(70000.0,width,ne_w), lockToGrid(20000.,height,ne)]
# s=[1.,0.,0.]

dom=Brick(l0=width, l1=width, l2=height, n0=ne_w, n1=ne_w, n2=ne)

print "fstart:",fstart
print "fend:",fend
print "slip ",s

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
     print i,p[i]
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
u=pde.getSolution(verbose=True)
g_s=grad(u)
sigma=mu*symmetric(g_s)+lmbd*trace(g_s)*kronecker(d)
saveVTK("dis.xml",disp=u,sigma=sigma,cfs=sigma[0,1]-0.4*sigma[0,0])
