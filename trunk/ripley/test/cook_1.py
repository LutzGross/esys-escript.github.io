#!/usr/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from esys.escript import *
from esys.finley import Rectangle
from esys.escript.linearPDEs import LinearSinglePDE


matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

NE=100
H=1.
L=1.


class PoissonTest(object):
     def __init__(self, domain, f0, r, xc):
          self.domain=domain
          DIM=self.domain.getDim()
	  self.pde=LinearSinglePDE(self.domain)
          x=self.domain.getX()
          self.pde.setValue(A=kronecker(domain), q=whereZero(x[DIM-1]-inf(x[DIM-1])))
          self.u_ref=self.getU(f0, r, xc)
          self.f0=f0
          self.r=r
          self.xc=xc
          self.omega=1.

     def getU(self, f0, r, xc):
          x=Function(self.domain).getX()
          self.pde.setValue(Y=exp(f0-length(x-xc)**2/r**2))
          u=self.pde.getSolution()
          print "evaluated for f0, r, xc, u=",f0, r, xc,inf(u),sup(u)
          return u

     def getCostFunction(self, f0, r, xc):
            u=interpolate(self.getU(f0, r, xc), FunctionOnBoundary(self.domain))
            u_ref=interpolate(self.u_ref, FunctionOnBoundary(self.domain))
            C=0.5*integrate(self.omega*length(u-u_ref)**2)
            return C

def plot(X, Y, Z, filename="test.png", title="None"):

    plt.figure()
    im = plt.imshow(Z, interpolation='bilinear', origin='lower',
                    cmap=cm.gray,
                    extent=(0,1.,0,1.))
    
    levels = np.linspace(np.amin(Z), np.amax(Z), 20)
    CS = plt.contour(X, Y, Z,levels,
                     origin='lower',
                     linewidths=2)

    #Thicken the zero contour.
    zc = CS.collections[6]
    plt.setp(zc, linewidth=4)

    plt.clabel(CS, levels[1::2],  # label every second level
               inline=1,
               fmt='%2.2e',
               fontsize=14)

    # make a colorbar for the contour lines
    #CB = plt.colorbar(CS, shrink=0.8, extend='both')

    plt.title(title)
    plt.xlabel('r')
    plt.ylabel('f_0')
    #plt.hot()  # Now change the colormap for the contour lines and colorbar
    plt.flag()

    # We can still add a colorbar for the image, too.
    # CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)

    # This makes the original colorbar look a bit out of place,
    # so let's improve its position.

    # l,b,w,h = plt.gca().get_position().bounds
    # ll,bb,ww,hh = CB.ax.get_position().bounds
    # CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])


    plt.savefig(filename, dpi=100)

def runTestOverGrid(p):
    delta=0.1
    # r = np.arange(max(delta,2./NE), 1.0+delta, delta)
    r = np.arange(0., 1.0+delta, delta)
    f0 = np.arange(0, 1.0+delta, delta)
    R, F0= np.meshgrid(r, f0)
    C=R.copy()
    for i in range(C.shape[0]):
       for j in range(C.shape[1]):
          # C[i,j]=p.getCostFunction(F0[i,j], R[i,j], p.xc)
          C[i,j]=p.getCostFunction(p.f0, p.r, (R[i,j],F0[i,j]))
    plot(R,F0,C, "test.png", "test")

p1=PoissonTest(Rectangle(NE,int(NE*H/L),l0=L,l1=H),  f0=0.3, r=0.1, xc=(0.5,H*0.5))
runTestOverGrid(p1)
