# $Id$
from esys.linearPDEs import LinearPDE
import numarray
class Helmholtz(LinearPDE):
   def setValue(self,kappa=0,omega=1,f=0,eta=0,g=0):
        self._setValue(A=kappa*numarray.identity(self.getDim()),D=omega,Y=f,d=eta,y=g)
