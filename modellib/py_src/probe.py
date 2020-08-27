
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript.modelframe import Model,ParameterSet
from esys.escript import Data
from esys.escript.util import *


class EvaluateExpression(ParameterSet):
       """
       Return the evaluation of an expression at current time t and 
       locations in the domain

       :warning: this class use python's eval function!!!!! 
                 Please use input.InterpolateOverBox is possible!!!!
       :note: Instance variable expression - expression or list of expressions defining 
                 expression value (in)
       :note: Instance variable out - current value of the expression (callable)
       """

       def __init__(self,**kwargs):
           """
           Set up parameters
           """
           super(EvaluateExpression, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 t=0., \
                                 expression="x[0]")

       def out(self):
         x=self.domain.getX()
         t=self.t
         if isinstance(self.expression,str):
           out=eval(self.expression)
         else:
           out=Data(0,(len(self.expression),),x.getFunctionSpace())
           for i in range(len(self.expression)): out[i]=eval(self.expression[i])
         return out

class Probe(Model):
       """
       Tests values against a expression which may depend on time and spatial 
       coordinates.
       
       It prints out the relative error in each time step and the maximum
       relative error over all time steps at the end.

       :warning: this class uses python's eval function!!!!!
       :ivar value: values to be tested (in)
       :ivar expression: expressions defining expression values to test against. If None only value is reported.  (in)
       :ivar line_tag: tag to be used when printing error.  (in)
       :ivar t: current time  (in)
       :ivar max_error: maximum error  (out)
       :ivar t_max: time of maximum error  (out)

       """

       def __init__(self,**kwargs):
           """
           Set up parameters
           """
           super(Probe, self).__init__(**kwargs)
           self.declareParameter(expression=None, \
                                 value=0., \
                                 t=0., \
                                 line_tag="PROBE")

       def doInitialization(self):
           """
           Initializes values
           """
           self.t_max=None
           self.max_error=0.

       def doStepPostprocessing(self,dt):
            t=self.t
            print("%s : time %e"%(self.line_tag,t))
            v=self.value
            x=v.getFunctionSpace().getX()
            if v.getRank()==0:
              if self.expression==None:
                 print("%s : true (min,max)= (%e,%e)"%(self.line_tag,inf(v),sub(v)))
              else:
                ref=eval(self.expression)
                err=Lsup(v-ref)/Lsup(ref)
                print("%s : true (min,max)= (%e,%e)"%(self.line_tag,inf(ref),sup(ref)))
                print("%s : (min,max,rel error)= (%e,%e,%e)"%(self.line_tag,inf(v),sup(v),err))
            else:
              err=0
              for i in range(v.getRank()):
                 vi=v[i]
                 if self.expression==None:
                    print("%s :   component %d: true (min,max)= (%e,%e)"%(self.line_tag,i,inf(vi)))
                 else:
                    refi=eval(self.expression[i])
                    erri=Lsup(vi-refi)/Lsup(refi)
                    print("%s :   component %d true (min,max)= (%e,%e)"%(self.line_tag,i,inf(refi),sup(refi)))
                    print("%s :   component %d (min,max,rel error)= (%e,%e,%e,%e,%e)"%(self.line_tag,i,inf(vi),sup(vi),erri))
                    err=max(err,erri)
              if not self.expression==None: print("%s :   maximum error %e",err)
            
            if not self.expression==None: 
               if err>self.max_error:
                   self.t_max=t
                   self.max_error=err

       def doFinalization(self):
          """
          Print out the maximum error.
          """
          if not self.t_max==None: print("%s : == maximum error %e at time %e == "%(self.line_tag,self.max_error,self.t_max))

# vim: expandtab shiftwidth=4:
