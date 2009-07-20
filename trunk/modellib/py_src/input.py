
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.modelframe import Model,ParameterSet
from esys.escript.linearPDEs import LinearPDE
from math import log

class Sequencer(Model):
    """
    Runs through time until t_end is reached. 

    @ivar t_end: model is terminated when t_end is passed, default 1 (in).
    @type t_end: C{float}
    @ivar dt_max: maximum time step size, default L{Model.UNDEF_DT} (in)
    @type dt_max: C{float}
    @ivar t: current time stamp (in/out). By default it is initialized with zero.
    @type t: C{float}

    """
    def __init__(self,**kwargs):
        """
        """
        super(Sequencer,self).__init__(**kwargs)
        self.declareParameter(t=0.,
                              t_end=1.,
                              dt_max=Model.UNDEF_DT)

    def doInitialization(self):
        """ 
        initialize time integration
        """
        self.__t_old = self.t

    def doStepPreprocessing(self, dt):
        self.t = self.__t_old+dt

    def doStepPostprocessing(self, dt):
        self.__t_old = self.t

    def finalize(self):
        """
        returns true when L{t} has reached L{t_end}
        """
        return self.t >= self.t_end

    def getSafeTimeStepSize(self, dt):
        """
        returns L{dt_max}
        """
        return self.dt_max

class GaussianProfile(ParameterSet):
    """
    Generates a Gaussian profile at center x_c, width width and height A 
    over a domain

    @ivar domain: domain
    @ivar x_c: center of the Gaussian profile (default [0.,0.,0.])
    @ivar A: (in) height of the profile. A maybe a vector. (default 1.)
    @ivar width: (in) width of the profile (default 0.1)
    @ivar r: (in) radius of the circle (default = 0)

    In the case that the spatial dimension is two, The third component of 
    x_c is dropped.
    """
    def __init__(self,**kwargs):
        super(GaussianProfile, self).__init__(**kwargs)
        self.declareParameter(domain=None, 
                              x_c=numpy.zeros([3]),
                              A=1.,
                              width=0.1,
                              r=0)

    def out(self):
        """
        Generate the Gaussian profile

        Link against this method to get the output of this model.
        """
        x = self.domain.getX()
        dim = self.domain.getDim()
        l = length(x-self.x_c[:dim])
        m = whereNegative(l-self.r)

        return (m+(1.-m)*exp(-log(2.)*(l/self.width)**2))*self.A

class InterpolateOverBox(ParameterSet):
    """
    Returns values at each time. The values are defined through given values 
    at time node. For two dimensional domains back values are ignored.

    @ivar domain: domain
    @ivar value_left_bottom_front: (in) value at left,bottom,front corner
    @ivar value_right_bottom_front: (in) value at right, bottom, front corner
    @ivar value_left_top_front: (in) value at left,top,front corner
    @ivar value_right_top_front: (in) value at right,top,front corner
    @ivar value_left_bottom_back: (in) value at  left,bottom,back corner
    @ivar value_right_bottom_back: (in) value at right,bottom,back corner
    @ivar value_left_top_back: (in) value at left,top,back  corner
    @ivar value_right_top_back: (in) value at right,top,back corner
    """

    def __init__(self, **kwargs):
        super(InterpolateOverBox, self).__init__(self)
        self.declareParameter(domain=None, 
                              value_left_bottom_front=0.,
                              value_right_bottom_front=0.,
                              value_left_top_front=0.,
                              value_right_top_front=0.,
                              value_left_bottom_back=0.,
                              value_right_bottom_back=0.,
                              value_left_top_back=0.,
                              value_right_top_back=0.)


    def out(self):
        """
        values at domain locations by bilinear interpolation of the given values.

        Link against this method to get the output of this model.
        """
        x = self.domain.getX()
        if self.domain.getDim() == 2:
            x0,x1=x[0],x[1]
            left_bottom_front0,right_top_back0=inf(x0),sup(x0)
            left_bottom_front1,right_top_back1=inf(x1),sup(x1)
            f_right = (x0 - left_bottom_front0)/(right_top_back0 -left_bottom_front0)
            f_left = 1. - f_right
            f_top = (x1 - left_bottom_front1)/(right_top_back1 - left_bottom_front1)
            f_bottom = 1. - f_top
            out = f_left * f_bottom * self.value_left_bottom_front \
                + f_right * f_bottom * self.value_right_bottom_front \
                + f_left * f_top * self.value_left_top_front \
                + f_right * f_top * self.value_right_top_front
        else:
            x0,x1,x2=x[0],x[1],x[2]
            left_bottom_front0,right_top_back0=inf(x0),sup(x0)
            left_bottom_front1,right_top_back1=inf(x1),sup(x1)
            left_bottom_front2,right_top_back2=inf(x2),sup(x2)
            f_right = (x0 - left_bottom_front0)/(right_top_back0 - left_bottom_front0)
            f_left = 1. - f_right
            f_top = (x1 - left_bottom_front1)/(right_top_back1 - left_bottom_front1)
            f_bottom = 1. - f_top
            f_back = (x2 - left_bottom_front1)/(right_top_back2 - left_bottom_front2)
            f_front = 1. - f_back
            out = f_left * f_bottom * f_front * self.value_left_bottom_front\
                + f_right * f_bottom * f_front * self.value_right_bottom_front\
                + f_left * f_top * f_front * self.value_left_top_front\
                + f_right * f_top * f_front * self.value_right_top_front\
                + f_left * f_bottom * f_back * self.value_left_bottom_back\
                + f_right * f_bottom * f_back * self.value_right_bottom_back\
                + f_left * f_top * f_back * self.value_left_top_back\
                + f_right * f_top * f_back * self.value_right_top_back
        return out


class InterpolatedTimeProfile(ParameterSet):
       """

       Returns values at each time. The values are defined through given 
       values at time node.
            
       value[i] defines the value at time nodes[i]. Between nodes linear 
       interpolation is used.

       For time t<nodes[0], value[0] is used and for t>nodes[l], values[l] 
       is used where l=len(nodes)-1.
 
       @ivar t: (in) current time
       @ivar node: (in) list of time nodes
       @ivar values: (in) list of values at time nodes
       """

       def __init__(self,**kwargs):
           super( InterpolatedTimeProfile, self).__init__(**kwargs)
           self.declareParameter(t=0., \
                                 nodes=[0.,1.],\
                                 values=[1.,1.])
       def out(self):
           """
           current value
  
           Link against this method to get the output of this model.
           """
           l = len(self.nodes) - 1
           t = self.t
           if t <= self.nodes[0]:
               return self.values[0]
           else:
               for i in range(1,l):
                  if t < self.nodes[i]:
                      m = (self.values[i-1] - self.values[i])/\
                            (self.nodes[i-1] - self.nodes[i])
                      return m*(t-self.nodes[i-1]) + self.values[i-1]
               return self.values[l]

class ScalarDistributionFromTags(ParameterSet):
    """
    creates a scalar distribution on a domain from tags, If tag_map is given
    the tags can be given a names and tag_map is used to map it into domain tags.
            
    @ivar domain: domain
    @type domain: L{esys.escript.Domain}
    @ivar default: default value 
    @ivar tag0: tag 0
    @type tag0: C{int}
    @ivar value0: value for tag 0
    @type value0: C{float}
    @ivar tag1: tag 1
    @type tag1: C{int}
    @ivar value1: value for tag 1
    @type value1: C{float}
    @ivar tag2: tag 2
    @type tag2: C{int}
    @ivar value2: value for tag 2
    @type value2: C{float}
    @ivar tag3: tag 3
    @type tag3: C{int}
    @ivar value3: value for tag 3
    @type value3: C{float}
    @ivar tag4: tag 4
    @type tag4: C{int}
    @ivar value4: value for tag 4
    @type value4: C{float}
    @ivar tag5: tag 5
    @type tag5: C{int}
    @ivar value5: value for tag 5
    @type value5: C{float}
    @ivar tag6: tag 6
    @type tag6: C{int}
    @ivar value6: value for tag 6
    @type value6: C{float}
    @ivar tag7: tag 7
    @type tag7: C{int}
    @ivar value7: value for tag 7
    @type value7: C{float}
    @ivar tag8: tag 8
    @type tag8: C{int}
    @ivar value8: value for tag 8
    @type value8: C{float}
    @ivar tag9: tag 9
    @type tag9: C{int}
    @ivar value9: value for tag 9
    @type value9: C{float}
    """
    def __init__(self,**kwargs):
        super(ScalarDistributionFromTags, self).__init__(**kwargs)
        self.declareParameter(domain=None,
                              default=0.,
                              tag0=None,
                              value0=0.,
                              tag1=None,
                              value1=0.,
                              tag2=None,
                              value2=0.,
                              tag3=None,
                              value3=0.,
                              tag4=None,
                              value4=0.,
                              tag5=None,
                              value5=0.,
                              tag6=None,
                              value6=0.,
                              tag7=None,
                              value7=0.,
                              tag8=None,
                              value8=0.,
                              tag9=None,
                              value9=0.)


    def out(self):
        """
        returns a L{esys.escript.Data} object
        Link against this method to get the output of this model.
        """
        d=Scalar(self.default,Function(self.domain))
        if not self.tag0 == None: d.setTaggedValue(self.tag0,self.value0)
        if not self.tag1 == None: d.setTaggedValue(self.tag1,self.value1)
        if not self.tag2 == None: d.setTaggedValue(self.tag2,self.value2)
        if not self.tag3 == None: d.setTaggedValue(self.tag3,self.value3)
        if not self.tag4 == None: d.setTaggedValue(self.tag4,self.value4)
        if not self.tag5 == None: d.setTaggedValue(self.tag5,self.value5)
        if not self.tag6 == None: d.setTaggedValue(self.tag6,self.value6)
        if not self.tag7 == None: d.setTaggedValue(self.tag7,self.value7)
        if not self.tag8 == None: d.setTaggedValue(self.tag8,self.value8)
        if not self.tag9 == None: d.setTaggedValue(self.tag9,self.value9)
        return d

class SmoothScalarDistributionFromTags(ParameterSet):
    """
    creates a smooth scalar distribution on a domain from region tags
            
    @ivar domain: domain
    @type domain: L{esys.escript.Domain}
    @ivar default: default value 
    @ivar tag0: tag 0
    @type tag0: C{int}
    @ivar value0: value for tag 0
    @type value0: C{float}
    @ivar tag1: tag 1
    @type tag1: C{int}
    @ivar value1: value for tag 1
    @type value1: C{float}
    @ivar tag2: tag 2
    @type tag2: C{int}
    @ivar value2: value for tag 2
    @type value2: C{float}
    @ivar tag3: tag 3
    @type tag3: C{int}
    @ivar value3: value for tag 3
    @type value3: C{float}
    @ivar tag4: tag 4
    @type tag4: C{int}
    @ivar value4: value for tag 4
    @type value4: C{float}
    @ivar tag5: tag 5
    @type tag5: C{int}
    @ivar value5: value for tag 5
    @type value5: C{float}
    @ivar tag6: tag 6
    @type tag6: C{int}
    @ivar value6: value for tag 6
    @type value6: C{float}
    @ivar tag7: tag 7
    @type tag7: C{int}
    @ivar value7: value for tag 7
    @type value7: C{float}
    @ivar tag8: tag 8
    @type tag8: C{int}
    @ivar value8: value for tag 8
    @type value8: C{float}
    @ivar tag9: tag 9
    @type tag9: C{int}
    @ivar value9: value for tag 9
    @type value9: C{float}
    """
    def __init__(self,**kwargs):
        super(SmoothScalarDistributionFromTags, self).__init__(**kwargs)
        self.declareParameter(domain=None,
                              default=0.,
                              tag0=None,
                              value0=0.,
                              tag1=None,
                              value1=0.,
                              tag2=None,
                              value2=0.,
                              tag3=None,
                              value3=0.,
                              tag4=None,
                              value4=0.,
                              tag5=None,
                              value5=0.,
                              tag6=None,
                              value6=0.,
                              tag7=None,
                              value7=0.,
                              tag8=None,
                              value8=0.,
                              tag9=None,
                              value9=0.)


    def __update(self,tag,tag_value,value):
        if self.__pde==None:
           self.__pde=LinearPDE(self.domain,numSolutions=1)
        mask=Scalar(0.,Function(self.domain))
        mask.setTaggedValue(tag,1.)
        self.__pde.setValue(Y=mask)
        mask=wherePositive(abs(self.__pde.getRightHandSide()))
        value*=(1.-mask)
        value+=tag_value*mask
        return value

    def out(self):
        """
        returns a L{esys.escript.Data} object
        Link against this method to get the output of this model.
        """
        d=Scalar(self.default,Solution(self.domain)) 
        self.__pde=None
        if not self.tag0 == None: d=self.__update(self.tag0,self.value0,d)
        if not self.tag1 == None: d=self.__update(self.tag1,self.value1,d)
        if not self.tag2 == None: d=self.__update(self.tag2,self.value2,d)
        if not self.tag3 == None: d=self.__update(self.tag3,self.value3,d)
        if not self.tag4 == None: d=self.__update(self.tag4,self.value4,d)
        if not self.tag5 == None: d=self.__update(self.tag5,self.value5,d)
        if not self.tag6 == None: d=self.__update(self.tag6,self.value6,d)
        if not self.tag7 == None: d=self.__update(self.tag7,self.value7,d)
        if not self.tag8 == None: d=self.__update(self.tag8,self.value8,d)
        if not self.tag9 == None: d=self.__update(self.tag9,self.value9,d)
        return d

class LinearCombination(ParameterSet):
    """
    Returns a linear combination of the f0*v0+f1*v1+f2*v2+f3*v3+f4*v4
            
    @ivar f0: numerical object or None, default=None (in)
    @ivar v0: numerical object or None, default=None (in)
    @ivar f1: numerical object or None, default=None (in)
    @ivar v1: numerical object or None, default=None (in)
    @ivar f2: numerical object or None, default=None (in)
    @ivar v2: numerical object or None, default=None (in)
    @ivar f3: numerical object or None, default=None (in)
    @ivar v3: numerical object or None, default=None (in)
    @ivar f4: numerical object or None, default=None (in)
    @ivar v4: numerical object or None, default=None (in)
    """
    def __init__(self,**kwargs):
        super(LinearCombination, self).__init__(**kwargs)
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
        """
        returns f0*v0+f1*v1+f2*v2+f3*v3+f4*v4.
        Link against this method to get the output of this model.
        """
        if not self.f0 == None and not self.v0 == None:
            fv0 = self.f0*self.v0
        else:
            fv0 = None

        if not self.f1 == None and not self.v1 == None:
            fv1 = self.f1*self.v1
        else:
            fv1 = None

        if not self.f2 == None and not self.v2 == None:
            fv2 = f2*v2
        else:
            fv2 = None

        if not self.f3 == None and not self.v3 == None:
            fv3 = self.f3*self.v3
        else:
            fv3 = None

        if not self.f4 == None and not self.v4 == None:
            fv4 = self.f4*self.v4
        else:
            fv4 = None

        if fv0 == None: 
             out = 0.
        else:
             out = fv0
        if not fv1 == None: 
            out += fv1
        if not fv2 == None: 
            out += fv2
        if not fv3 == None: 
            out += fv3
        return out

class MergeConstraints(ParameterSet):
    """
    Returns a linear combination of the f0*v0+f1*v1+f2*v2+f3*v3+f4*v4
    """
    def __init__(self,**kwargs):
        super(MergeConstraints, self).__init__(**kwargs)
        self.declareParameter(location_of_constraint0=None, \
                              value_of_constraint0=None, \
                              location_of_constraint1=None, \
                              value_of_constraint1=None, \
                              location_of_constraint2=None, \
                              value_of_constraint2=None, \
                              location_of_constraint3=None, \
                              value_of_constraint3=None, \
                              location_of_constraint4=None, \
                              value_of_constraint4=None)
    def location_of_constraint(self):
          """
          return the values used to constrain a solution

          @return: the mask marking the locations of the constraints
          @rtype: L{escript.Scalar}
          """
          out_loc=0
          if not self.location_of_constraint0 == None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint0))
          if not self.location_of_constraint1 == None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint1))
          if not self.location_of_constraint2 == None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint2))
          if not self.location_of_constraint3 == None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint3))
          return out_loc

    def value_of_constraint(self):
          """
          return the values used to constrain a solution

          @return: values to be used at the locations of the constraints. If
                  C{value} is not given C{None} is rerturned.
          @rtype: L{escript.Scalar}
          """
          out_loc=0
          out=0
          if not self.location_of_constraint0 == None:
               tmp=wherePositive(self.location_of_constraint0)
               out=out*(1.-tmp)+self.value_of_constraint0*tmp
               out_loc=wherePositive(out_loc+tmp)
          if not self.location_of_constraint1 == None:
               tmp=wherePositive(self.location_of_constraint1)
               out=out*(1.-tmp)+self.value_of_constraint1*tmp
               out_loc=wherePositive(out_loc+tmp)
          if not self.location_of_constraint2 == None:
               tmp=wherePositive(self.location_of_constraint2)
               out=out*(1.-tmp)+self.value_of_constraint2*tmp
               out_loc=wherePositive(out_loc+tmp)
          if not self.location_of_constraint3 == None:
               tmp=wherePositive(self.location_of_constraint3)
               out=out*(1.-tmp)+self.value_of_constraint3*tmp
               out_loc=wherePositive(out_loc+tmp)
          return out
# vim: expandtab shiftwidth=4:
