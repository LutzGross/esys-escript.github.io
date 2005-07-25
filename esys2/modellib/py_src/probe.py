# $Id$


from esys.modelframe import Model
from esys.escript import inf,sup,Lsup


class Probe(Model):
       """@brief tests values against a reference which may depend on time and spatial coordinates
                 it prints out the relative error in each time step and the maximum relative error over
                 all time steps at the end

           WARNING: this class use python's eval function!!!!!

          @param value (in) - values to be tested
          @param reference (in) - expressions defining reference values to test against. If None only value is reported.
          @param line_tag (in) - tag to be used when printing error
          @param max_error (out) - maximum error
          @param t_max (out) - time of maximum error

        """

       def __init__(self,debug=False):
           """set up parameters"""
           Model.__init__(self,debug=debug)
           self.declareParameter(reference=None, \
                                 value=0., \
                                 line_tag="PROBE")

       def doInitialization(self,t):
           """initializes values"""
           self.__tn=t
           self.t_max=None
           self.max_err=0.

       def doStep(self,dt):
            t=self.__tn+dt
            print "%s : time %d"%t
            v=self.value
            x=v.getFunctionSpace().getX()
            if v.getRank()==0:
              if self.reference==None:
                 print "%s :   (min,max)= (%e,%e)"%(self.line_tag,inf(v),max(v))
              else:
                ref=eval(self.reference)
                err=Lsup(v-ref)/Lsup(ref)
                print "%s :   (min,max,rel error)= (%e,%e,%e)"%(self.line_tag,inf(v),max(v),err)
            else:
              err=0
              for i in range(v.getRank()):
                 vi=v[i]
                 if self.reference==None:
                    print "%s :   component %d (min,max)= (%e,%e)"%(self.line_tag,i,inf(vi),max(vi))
                 else:
                    refi=eval(self.reference[i])
                    erri=Lsup(vi-refi)/Lsup(refi)
                    print "%s :   component %d (min,max,rel error)= (%e,%e,%e)"%(self.line_tag,i,inf(vi),max(vi),erri)
                    err=max(err,erri)
              if not self.reference==None: print "%s :   maximum error %e",err
            self.__tn=t
            
            if not self.reference==None: 
               if err>self.max_err:
                   self.t_max=self.__tn
                   self.max_error=err/max_ref

       def finalize(self):
          """print out the maximum error"""
          if not self.t_max==None: print "%s :   component %d (min,max)= (%e,%e)"%(self.line_tag,self.t_max,self.max_err)
