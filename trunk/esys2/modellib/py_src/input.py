# $Id$
from escript.escript import *
from escript.modelframe import Model
from math import log

class Sequencer(Model):
    """@brief runs through time until t_end is reached.

           @param t_end (in) -  model is terminate when t_end is passed (default Model.UNDEF_DT)
           @param dt_max (in) -  maximum time step size (default Model.UNDEF_DT)
           @param t (out) -  current time

    """
    def __init__(self,debug=False):
        Model.__init__(self,debug=debug)
        self.declareParameter(t=0., \
                              t_end=Model.UNDEF_DT,  \
                              dt_max=Model.UNDEF_DT)

    def doInitialization(self):
      self.__t_old=self.t

    def doStepPreprocessing(self,dt):
      self.t=self.__t_old+dt

    def doStepPostprocessing(self,dt):
      self.__t_old=self.t

    def finalize(self):
      """true when t has reached t_end"""
      return self.t>=self.t_end

    def getSafeTimeStepSize(self,dt):
      """returns dt_max"""
      return self.dt_max

class GausseanProfile(Model):
    """@brief generates a gaussean profile at center x_c, width width and height A over a domain

           @param domain (in) - domain
           @param x_c (in)  - center of the Gaussean profile (default [0.,0.,0.])
           @param A (in)  - height of the profile. A maybe a vector. (default 1.)
           @param width (in) - width of the profile (default 0.1)
           @param r (in) -  radius of the circle (default = 0)
           @param out (callable) - profile


       In the case that the spatial dimension is two, The third component of x_c is dropped
    """
    def __init__(self,debug=False):
        Model.__init__(self,debug=debug)
        self.declareParameter(domain=None, 
                              x_c=numarray.zeros([3]),
                              A=1.,
                              width=0.1,
                              r=0)

    def out(self):
        x=self.domain.getX()
        dim=self.domain.getDim()
        l=length(x-self.x_c[:dim])
        m=(l-self.r).whereNegative()
        return (m+(1.-m)*exp(-log(2.)*(l/self.width)**2))*self.A

class InterpolateOverBox(Model):
    """
           @brief returns values at each time. The values are defined through given values at time node.

           @param domain (in) - domain
           @param left_bottom_front (in) - coordinates of left,bottom,front corner of the box
           @param right_top_back (in) - coordinates of the right, top, back corner of the box
           @param value_left_bottom_front (in) - value at left,bottom,front corner
           @param value_right_bottom_front (in) - value at right, bottom, front corner
           @param value_left_top_front (in) - value at left,top,front corner
           @param value_right_top_front (in) - value at right,top,front corner
           @param value_left_bottom_back (in) - value at  left,bottom,back corner
           @param value_right_bottom_back (in) - value at right,bottom,back corner
           @param value_left_top_back (in) - value at left,top,back  corner
           @param value_right_top_back (in) - value at right,top,back corner
           @param out (callable) - values at doamin locations by bilinear interpolation. for two dimensional domains back values are ignored.

    """

    def __init__(self,debug=False):
        Model.__init__(self,debug=debug)
        self.declareParameter(domain=None, 
                              left_bottom_front=[0.,0.,0.],
                              right_top_back=[1.,1.,1.],
                              value_left_bottom_front=0.,
                              value_right_bottom_front=0.,
                              value_left_top_front=0.,
                              value_right_top_front=0.,
                              value_left_bottom_back=0.,
                              value_right_bottom_back=0.,
                              value_left_top_back=0.,
                              value_right_top_back=0.)


    def out(self):
        x=self.domain.getX()
        if self.domain.getDim()==2:
         f_right=(x[0]-self.left_bottom_front[0])/(self.right_top_back[0]-self.left_bottom_front[0])
         f_left=1.-f_right
         f_top=(x[1]-self.left_bottom_front[1])/(self.right_top_back[1]-self.left_bottom_front[1])
         f_bottom=1.-f_top
         out=self.value_left_bottom_front * f_left * f_bottom \
            +self.value_right_bottom_front* f_right * f_bottom \
            +self.value_left_top_front    * f_left * f_top \
            +self.value_right_top_front   * f_right * f_top 

        else:
         f_right=(x[0]-self.left_bottom_front[0])/(self.right_top_back[0]-self.left_bottom_front[0])
         f_left=1.-f_right
         f_top=(x[1]-self.left_bottom_front[1])/(self.right_top_back[1]-self.left_bottom_front[1])
         f_bottom=1.-f_top
         f_back=(x[2]-self.left_bottom_front[1])/(self.right_top_back[2]-self.left_bottom_front[2])
         f_front=1.-f_back
         out=self.value_left_bottom_front * f_left * f_bottom * f_front \
            +self.value_right_bottom_front* f_right * f_bottom * f_front \
            +self.value_left_top_front    * f_left * f_top * f_front \
            +self.value_right_top_front   * f_right * f_top * f_front \
            +self.value_left_bottom_back  * f_left * f_bottom * f_back \
            +self.value_right_bottom_back * f_right * f_bottom * f_back \
            +self.value_left_top_back     * f_left * f_top * f_back \
            +self.value_right_top_back    * f_right * f_top * f_back
        return out


class InterpolatedTimeProfile(Model):
       """@brief returns values at each time. The values are defined through given values at time node.
            
         value[i] defines the value at time nodes[i]. Between nodes linear interpolation is used.
         for time t<nodes[0] value[0] is used and for t>nodes[l] values[l] is used where l=len(nodes)-1.
 

          @param t (in) - current time
          @param node (in) - list of time nodes
          @param values (in) - list of values at time nodes
          @param out (callable) - current value 


        """

       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(t=0., \
                                 nodes=[0.,1.],\
                                 values=[1.,1.])
       def out(self):
           l=len(self.nodes)-1
           t=self.t
           if t<=self.nodes[0]:
                return self.values[0]
           else:
               for i in range(1,l):
                  if t<self.nodes[i]:
                    m=(self.values[i-1]-self.values[i])/(self.nodes[i-1]-self.nodes[i])
                    return m*(t-self.nodes[i-1])+self.values[i-1]
               return self.values[l]

class LinearCombination(Model):
       """@brief returns a linear combination of the f0*v0+f1*v1+f2*v2+f3*v3+f4*v4
            
         @param f0 (in) numerical object or None (default: None) 
         @param v0 (in) numerical object or None (default: None) 
         @param f1 (in) numerical object or None (default: None) 
         @param v1 (in) numerical object or None (default: None) 
         @param f2 (in) numerical object or None (default: None) 
         @param v2 (in) numerical object or None (default: None) 
         @param f3 (in) numerical object or None (default: None) 
         @param v3 (in) numerical object or None (default: None) 
         @param f4 (in) numerical object or None (default: None) 
         @param v4 (in) numerical object or None (default: None)
         @param out (callable) - current value 


        """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(f0=None, \
                                 v0=None, \
                                 f1=None, \
                                 v1=None, \
                                 f2=None, \
                                 v2=None, \
                                 f3=None, \
                                 v3=None, \
                                 f4=None, \
                                 v4=None)
       def out(self):
           if not self.f0==None and not self.v0==None:
                fv0=self.f0*self.v0
           else:
                fv0=None

           if not self.f1==None and not self.v1==None:
                fv1=self.f1*self.v1
           else:
                fv1=None

           if not self.f2==None and not self.v2==None:
                fv2=f2*v2
           else:
                fv2=None

           if not self.f3==None and not self.v3==None:
                fv3=self.f3*self.v3
           else:
                fv3=None

           if not self.f4==None and not self.v4==None:
                fv4=self.f4*self.v4
           else:
                fv4=None

           if fv0==None: 
              out=0.
           else:
              out=fv0
           if not fv1==None: out+=fv1
           if not fv2==None: out+=fv2
           if not fv3==None: out+=fv3
           return out
