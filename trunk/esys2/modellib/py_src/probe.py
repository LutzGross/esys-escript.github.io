# $Id$

from escript.modelframe import Model,ParameterSet
from escript.escript import Data
from escript.util import *


class EvaluateExpression(ParameterSet):
       """@brief return the evaluation of an expression at current time t and locations in the domain

           WARNING: this class use python's eval function!!!!! Please use input.InterpolateOverBox is possible!!!!

          @param expression (in) - expression or list of expressions defining expression value
          @param out (out) - current value of the expression

        """

       def __init__(self,debug=False):
           """set up parameters"""
           ParameterSet.__init__(self,debug=debug)
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
       """@brief tests values against a expression which may depend on time and spatial coordinates
                 it prints out the relative error in each time step and the maximum relative error over
                 all time steps at the end

           WARNING: this class use python's eval function!!!!!

          @param value (in) - values to be tested
          @param expression (in) - expressions defining expression values to test against. If None only value is reported.
          @param line_tag (in) - tag to be used when printing error
          @param t (in) - current time
          @param max_error (out) - maximum error
          @param t_max (out) - time of maximum error

        """

       def __init__(self,debug=False):
           """set up parameters"""
           Model.__init__(self,debug=debug)
           self.declareParameter(expression=None, \
                                 value=0., \
                                 t=0., \
                                 line_tag="PROBE")

       def doInitialization(self):
           """initializes values"""
           self.t_max=None
           self.max_error=0.

       def doStepPostprocessing(self,dt):
            t=self.t
            print "%s : time %e"%(self.line_tag,t)
            v=self.value
            x=v.getFunctionSpace().getX()
            if v.getRank()==0:
              if self.expression==None:
                 print "%s : true (min,max)= (%e,%e)"%(self.line_tag,inf(v),sub(v))
              else:
                ref=eval(self.expression)
                err=Lsup(v-ref)/Lsup(ref)
                print "%s : true (min,max)= (%e,%e)"%(self.line_tag,inf(ref),sup(ref))
                print "%s : (min,max,rel error)= (%e,%e,%e)"%(self.line_tag,inf(v),sup(v),err)
            else:
              err=0
              for i in range(v.getRank()):
                 vi=v[i]
                 if self.expression==None:
                    print "%s :   component %d: true (min,max)= (%e,%e)"%(self.line_tag,i,inf(vi))
                 else:
                    refi=eval(self.expression[i])
                    erri=Lsup(vi-refi)/Lsup(refi)
                    print "%s :   component %d true (min,max)= (%e,%e)"%(self.line_tag,i,inf(refi),sup(refi))
                    print "%s :   component %d (min,max,rel error)= (%e,%e,%e,%e,%e)"%(self.line_tag,i,inf(vi),sup(vi),erri)
                    err=max(err,erri)
              if not self.expression==None: print "%s :   maximum error %e",err
            
            if not self.expression==None: 
               if err>self.max_error:
                   self.t_max=t
                   self.max_error=err

       def doFinalization(self):
          """print out the maximum error"""
          if not self.t_max==None: print "%s : == maximum error %e at time %e == "%(self.line_tag,self.max_error,self.t_max)
