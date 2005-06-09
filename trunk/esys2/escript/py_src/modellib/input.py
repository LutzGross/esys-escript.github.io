# $Id$


from esys.modelframe import Model
from esys.escript import *
from math import log


class GausseanProfile(Model):
    """@brief generates a gaussean profile at center x_c, width width and height A over a domain

           @param domain (in) - domain
           @param x_c (in)  - center of the Gaussean profile (default [0.,0.,0.])
           @param A (in)  - height of the profile. A maybe a vector. (default 1.)
           @param width (in) - width of the profile (default 0.1)
           @param r (in) -  radius of the circle (default = 0)
           @param out (out) - profile


       In the case that the spatial dimension is two, The third component of x_c is dropped
    """
    def __init__(self,debug=False):
        Model.__init__(self,debug=debug)
        self.declareParameter(domain=None, x_c=numarray.zeros([3]),A=1.,width=0.1,r=0)

    def out(self):
        x=self.domain.getX()
        dim=self.domain.getDim()
        l=length(x-self.x_c[:dim])
        m=(l-self.r).whereNegative()
        return (m+(1.-m)*exp(-log(2.)*(l/self.width)**2))*self.A

class InterpolatedTimeProfile(Model):
       """ """

       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(t=[0.,1.],\
                                 values=[1.,1.],\
                                 out=0.)
       def doInitialization(self,t):
           self.__tn=t

       def doStep(self,dt):
            t=self.__tn+dt
            if t<=self.t[0]:
                self.out=self.values[0]
            else:
               for i in range(1,len(self.t)):
                  if t<self.t[i]:
                    m=(self.values[i-1]-self.values[i])/(self.t[i-1]-self.t[i])
                    self.out=m*(t-self.t[i-1])+self.values[i-1]
               self.out=self.t[len(self.t)-1]
            self.__tn+=dt
