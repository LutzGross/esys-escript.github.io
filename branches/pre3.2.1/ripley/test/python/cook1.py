import numpy as np
import matplotlib.pyplot as plt
import random

from esys.escript import * 
from esys.escript.linearPDEs import LinearPDE
from esys.escript.pdetools import Locator
from esys.finley import Rectangle,Brick 

import scipy.optimize as so

NE=200
NE=100
H=1.
padding=0.0
numContours=25
numPlotsGrid=25

class PoissonTest(object):
      def __init__(self, domain, a, b, xc):
          self.domain=domain
          self.a=a
          self.b=b
          self.xc=xc
          self.pde=LinearPDE(domain)
          self.pde.setSymmetryOn()
          x=self.domain.getX()
          DIM=domain.getDim()
          self.pde.setValue(A=kronecker(DIM), q=whereZero(x[DIM-1]-inf(x[DIM-1])))
          self.u_ref=self.getU(a, b, xc)
          self.omega=1.

      def getU(self, a, b, xc):
          x=Function(self.domain).getX()
          self.pde.setValue(Y=exp(a-(b*length(x-xc))**2))
          u=self.pde.getSolution()
          print "PoissonTest.getU: a,b,xc, u =",a,b,xc,inf(u),sup(u)
          return u

      def getCostFunction(self, a, b, xc):
          u=interpolate(self.getU(a, b, xc), FunctionOnBoundary(self.domain))
          C=0.5*integrate( self.omega * (u-self.u_ref)**2)
          return C

class PoissonTestTopData(PoissonTest):
       def __init__(self, domain, a, b, xc):
           PoissonTest.__init__(self, domain, a, b, xc)
           x=FunctionOnBoundary(self.domain).getX()
           DIM=self.domain.getDim()
           self.omega=whereZero(x[DIM-1]-sup(x[DIM-1]))
class PoissonTestDiscreteData(PoissonTest):
       def __init__(self, N, sigma, domain, a, b, xc):
             PoissonTest.__init__(self, domain, a, b, xc)
             self.N=N
             self.sigma=sigma
             x=domain.getX()
             x_min=inf(x[0])
             x_max=sup(x[0])
             z_max=sup(x[self.domain.getDim()-1])
             dx=(x_max-x_min)/N
             self.loc=Locator(self.domain, [ (dx/2+dx*i, z_max) for i in xrange(N)  ])
             uu=self.loc(self.u_ref)
             s=self.sigma*(max(uu)+min(uu))/2.
             print self.loc.getX()
             self.data = numpy.array([    random.gauss(d, s) for d in self.loc(self.u_ref) ])
             print self.data
       def getCostFunction(self, a, b, xc):
          u=numpy.array(self.loc(self.getU(a, b, xc)))
          C=(0.5/self.N)*numpy.sum((u-self.data)**2)
          return C
             
         

def plotContour(p, filename, title):
    x0=np.linspace(-padding, 1+padding, numPlotsGrid)
    X0, Y0 = np.meshgrid(x0, H*x0)

    plt.figure()
    C=X0.copy()
    for i in range(C.shape[0]):
       for j in range(C.shape[1]):
          C[i,j]=p.getCostFunction(p.a,p.b, (X0[i,j], Y0[i,j]))

    C_max=np.amax(C)
    C_min=np.amin(C)
    levels = np.linspace(C_min, C_max, numContours)
    cset1 = plt.contourf(X0, Y0, C, levels, cmap=plt.cm.get_cmap('jet', len(levels)-1))
    cset2 = plt.contour(X0, Y0, C, cset1.levels, colors = 'k', hold='on')
    for c in cset2.collections: c.set_linestyle('solid')
    plt.title('Cost function (%s)'%title)
    plt.xlabel("xc_0")
    plt.ylabel("xc_1")
    plt.colorbar(cset1)
    plt.savefig(filename)
          


plotContour(PoissonTestDiscreteData(5, 0.10, Rectangle(NE, int(NE*H),l1=H), 0,15,(0.5,0.5)), "cost_passion_discrete_w_cc.png","b=10, xc=(0.5,0.2)")

# plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.5,0.2)), "cost_passion_top_w_cb.png","b=10, xc=(0.5,0.2)")
# plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.5,0.5)), "cost_passion_top_w_cc.png","b=10, xc=(0.5,0.5)")
# plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.5,0.8)), "cost_passion_top_w_ct.png","b=10, xc=(0.5,0.8)")
# plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.2,0.2)), "cost_passion_top_w_lb.png","b=10, xc=(0.2,0.2)")
# plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.2,0.5)), "cost_passion_top_w_lc.png","b=10, xc=(0.2,0.5)")
# plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.2,0.8)), "cost_passion_top_w_lt.png","b=10, xc=(0.2,0.8)")
# plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.5,0.2)), "cost_passion_all_w_cb.png","b=10, xc=(0.5,0.2)")
# plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.5,0.5)), "cost_passion_all_w_cc.png","b=10, xc=(0.5,0.5)")
# plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.5,0.8)), "cost_passion_all_w_ct.png","b=10, xc=(0.5,0.8)")
# plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.2,0.2)), "cost_passion_all_w_lb.png","b=10, xc=(0.2,0.2)")
# plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.2,0.5)), "cost_passion_all_w_lc.png","b=10, xc=(0.2,0.5)")
# plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,10,(0.2,0.8)), "cost_passion_all_w_lt.png","b=10, xc=(0.2,0.8)")

#plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.5,0.2)), "cost_passion_top_n_cb.png","b=100, xc=(0.5,0.2)")
#plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.5,0.5)), "cost_passion_top_n_cc.png","b=100, xc=(0.5,0.5)")
#plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.5,0.8)), "cost_passion_top_n_ct.png","b=100, xc=(0.5,0.8)")
#plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.2,0.2)), "cost_passion_top_n_lb.png","b=100, xc=(0.2,0.2)")
#plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.2,0.5)), "cost_passion_top_n_lc.png","b=100, xc=(0.2,0.5)")
#plotContour(PoissonTestTopData(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.2,0.8)), "cost_passion_top_n_lt.png","b=100, xc=(0.2,0.8)")
#plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.5,0.2)), "cost_passion_all_n_cb.png","b=100, xc=(0.5,0.2)")
#plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.5,0.5)), "cost_passion_all_n_cc.png","b=100, xc=(0.5,0.5)")
#plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.5,0.8)), "cost_passion_all_n_ct.png","b=100, xc=(0.5,0.8)")
#plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.2,0.2)), "cost_passion_all_n_lb.png","b=100, xc=(0.2,0.2)")
#plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.2,0.5)), "cost_passion_all_n_lc.png","b=100, xc=(0.2,0.5)")
#plotContour(PoissonTest(Rectangle(NE, int(NE*H),l1=H), 0,100,(0.2,0.8)), "cost_passion_all_n_lt.png","b=100, xc=(0.2,0.8)")





