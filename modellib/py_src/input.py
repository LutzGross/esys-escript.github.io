
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import length, wherePositive, whereNegative, exp, inf, sup
from esys.escript.modelframe import Model,ParameterSet
from esys.escript.linearPDEs import LinearPDE
from math import log
import numpy

class Sequencer(Model):
    """
    Runs through time until t_end is reached. 

    :ivar t_end: model is terminated when t_end is passed, default 1 (in).
    :type t_end: ``float``
    :ivar dt_max: maximum time step size, default `Model.UNDEF_DT` (in)
    :type dt_max: ``float``
    :ivar t: current time stamp (in/out). By default it is initialized with zero.
    :type t: ``float``

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
        returns true when `t` has reached `t_end`
        """
        return self.t >= self.t_end

    def getSafeTimeStepSize(self, dt):
        """
        returns `dt_max`
        """
        return self.dt_max

class GaussianProfile(ParameterSet):
    """
    Generates a Gaussian profile at center x_c, width width and height A 
    over a domain

    :note: Instance variable domain - domain
    :note: Instance variable x_c - center of the Gaussian profile (default [0.,0.,0.])
    :note: Instance variable A - (in) height of the profile. A maybe a vector. (default 1.)
    :note: Instance variable width - (in) width of the profile (default 0.1)
    :note: Instance variable r - (in) radius of the circle (default = 0)

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

    :note: Instance variable domain - domain
    :note: Instance variable value_left_bottom_front - (in) value at left,bottom,front corner
    :note: Instance variable value_right_bottom_front - (in) value at right, bottom, front corner
    :note: Instance variable value_left_top_front - (in) value at left,top,front corner
    :note: Instance variable value_right_top_front - (in) value at right,top,front corner
    :note: Instance variable value_left_bottom_back - (in) value at  left,bottom,back corner
    :note: Instance variable value_right_bottom_back - (in) value at right,bottom,back corner
    :note: Instance variable value_left_top_back - (in) value at left,top,back  corner
    :note: Instance variable value_right_top_back - (in) value at right,top,back corner
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
 
       :note: Instance variable t - (in) current time
       :note: Instance variable node - (in) list of time nodes
       :note: Instance variable values - (in) list of values at time nodes
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
            
    :ivar domain: domain
    :type domain: `esys.escript.Domain`
    :ivar default: default value 
    :ivar tag0: tag 0
    :type tag0: ``int``
    :ivar value0: value for tag 0
    :type value0: ``float``
    :ivar tag1: tag 1
    :type tag1: ``int``
    :ivar value1: value for tag 1
    :type value1: ``float``
    :ivar tag2: tag 2
    :type tag2: ``int``
    :ivar value2: value for tag 2
    :type value2: ``float``
    :ivar tag3: tag 3
    :type tag3: ``int``
    :ivar value3: value for tag 3
    :type value3: ``float``
    :ivar tag4: tag 4
    :type tag4: ``int``
    :ivar value4: value for tag 4
    :type value4: ``float``
    :ivar tag5: tag 5
    :type tag5: ``int``
    :ivar value5: value for tag 5
    :type value5: ``float``
    :ivar tag6: tag 6
    :type tag6: ``int``
    :ivar value6: value for tag 6
    :type value6: ``float``
    :ivar tag7: tag 7
    :type tag7: ``int``
    :ivar value7: value for tag 7
    :type value7: ``float``
    :ivar tag8: tag 8
    :type tag8: ``int``
    :ivar value8: value for tag 8
    :type value8: ``float``
    :ivar tag9: tag 9
    :type tag9: ``int``
    :ivar value9: value for tag 9
    :type value9: ``float``
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
        returns a `esys.escript.Data` object
        Link against this method to get the output of this model.
        """
        d=Scalar(self.default,Function(self.domain))
        if not self.tag0 is None: d.setTaggedValue(self.tag0,self.value0)
        if not self.tag1 is None: d.setTaggedValue(self.tag1,self.value1)
        if not self.tag2 is None: d.setTaggedValue(self.tag2,self.value2)
        if not self.tag3 is None: d.setTaggedValue(self.tag3,self.value3)
        if not self.tag4 is None: d.setTaggedValue(self.tag4,self.value4)
        if not self.tag5 is None: d.setTaggedValue(self.tag5,self.value5)
        if not self.tag6 is None: d.setTaggedValue(self.tag6,self.value6)
        if not self.tag7 is None: d.setTaggedValue(self.tag7,self.value7)
        if not self.tag8 is None: d.setTaggedValue(self.tag8,self.value8)
        if not self.tag9 is None: d.setTaggedValue(self.tag9,self.value9)
        return d

class SmoothScalarDistributionFromTags(ParameterSet):
    """
    creates a smooth scalar distribution on a domain from region tags
            
    :ivar domain: domain
    :type domain: `esys.escript.Domain`
    :ivar default: default value 
    :ivar tag0: tag 0
    :type tag0: ``int``
    :ivar value0: value for tag 0
    :type value0: ``float``
    :ivar tag1: tag 1
    :type tag1: ``int``
    :ivar value1: value for tag 1
    :type value1: ``float``
    :ivar tag2: tag 2
    :type tag2: ``int``
    :ivar value2: value for tag 2
    :type value2: ``float``
    :ivar tag3: tag 3
    :type tag3: ``int``
    :ivar value3: value for tag 3
    :type value3: ``float``
    :ivar tag4: tag 4
    :type tag4: ``int``
    :ivar value4: value for tag 4
    :type value4: ``float``
    :ivar tag5: tag 5
    :type tag5: ``int``
    :ivar value5: value for tag 5
    :type value5: ``float``
    :ivar tag6: tag 6
    :type tag6: ``int``
    :ivar value6: value for tag 6
    :type value6: ``float``
    :ivar tag7: tag 7
    :type tag7: ``int``
    :ivar value7: value for tag 7
    :type value7: ``float``
    :ivar tag8: tag 8
    :type tag8: ``int``
    :ivar value8: value for tag 8
    :type value8: ``float``
    :ivar tag9: tag 9
    :type tag9: ``int``
    :ivar value9: value for tag 9
    :type value9: ``float``
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
        returns a `esys.escript.Data` object
        Link against this method to get the output of this model.
        """
        d=Scalar(self.default,Solution(self.domain)) 
        self.__pde=None
        if not self.tag0 is None: d=self.__update(self.tag0,self.value0,d)
        if not self.tag1 is None: d=self.__update(self.tag1,self.value1,d)
        if not self.tag2 is None: d=self.__update(self.tag2,self.value2,d)
        if not self.tag3 is None: d=self.__update(self.tag3,self.value3,d)
        if not self.tag4 is None: d=self.__update(self.tag4,self.value4,d)
        if not self.tag5 is None: d=self.__update(self.tag5,self.value5,d)
        if not self.tag6 is None: d=self.__update(self.tag6,self.value6,d)
        if not self.tag7 is None: d=self.__update(self.tag7,self.value7,d)
        if not self.tag8 is None: d=self.__update(self.tag8,self.value8,d)
        if not self.tag9 is None: d=self.__update(self.tag9,self.value9,d)
        return d

class LinearCombination(ParameterSet):
    """
    Returns a linear combination of the f0*v0+f1*v1+f2*v2+f3*v3+f4*v4
            
    :ivar f0: numerical object or None, default=None (in)
    :ivar v0: numerical object or None, default=None (in)
    :ivar f1: numerical object or None, default=None (in)
    :ivar v1: numerical object or None, default=None (in)
    :ivar f2: numerical object or None, default=None (in)
    :ivar v2: numerical object or None, default=None (in)
    :ivar f3: numerical object or None, default=None (in)
    :ivar v3: numerical object or None, default=None (in)
    :ivar f4: numerical object or None, default=None (in)
    :ivar v4: numerical object or None, default=None (in)
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
        if not self.f0 is None and not self.v0 is None:
            fv0 = self.f0*self.v0
        else:
            fv0 = None

        if not self.f1 is None and not self.v1 is None:
            fv1 = self.f1*self.v1
        else:
            fv1 = None

        if not self.f2 is None and not self.v2 is None:
            fv2 = f2*v2
        else:
            fv2 = None

        if not self.f3 is None and not self.v3 is None:
            fv3 = self.f3*self.v3
        else:
            fv3 = None

        if not self.f4 is None and not self.v4 is None:
            fv4 = self.f4*self.v4
        else:
            fv4 = None

        if fv0 is None: 
             out = 0.
        else:
             out = fv0
        if not fv1 is None: 
            out += fv1
        if not fv2 is None: 
            out += fv2
        if not fv3 is None: 
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

          :return: the mask marking the locations of the constraints
          :rtype: `escript.Scalar`
          """
          out_loc=0
          if not self.location_of_constraint0 is None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint0))
          if not self.location_of_constraint1 is None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint1))
          if not self.location_of_constraint2 is None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint2))
          if not self.location_of_constraint3 is None:
               out_loc=wherePositive(out_loc+wherePositive(self.location_of_constraint3))
          return out_loc

    def value_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: values to be used at the locations of the constraints. If
                  ``value`` is not given ``None`` is rerturned.
          :rtype: `escript.Scalar`
          """
          out_loc=0
          out=0
          if not self.location_of_constraint0 is None:
               tmp=wherePositive(self.location_of_constraint0)
               out=out*(1.-tmp)+self.value_of_constraint0*tmp
               out_loc=wherePositive(out_loc+tmp)
          if not self.location_of_constraint1 is None:
               tmp=wherePositive(self.location_of_constraint1)
               out=out*(1.-tmp)+self.value_of_constraint1*tmp
               out_loc=wherePositive(out_loc+tmp)
          if not self.location_of_constraint2 is None:
               tmp=wherePositive(self.location_of_constraint2)
               out=out*(1.-tmp)+self.value_of_constraint2*tmp
               out_loc=wherePositive(out_loc+tmp)
          if not self.location_of_constraint3 is None:
               tmp=wherePositive(self.location_of_constraint3)
               out=out*(1.-tmp)+self.value_of_constraint3*tmp
               out_loc=wherePositive(out_loc+tmp)
          return out
# vim: expandtab shiftwidth=4:
