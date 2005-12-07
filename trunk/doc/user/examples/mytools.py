# $Id$
from esys.escript.linearPDEs import LinearPDE
from esys.escript import kronecker
import numarray
class Helmholtz(LinearPDE):
   def setValue(self,kappa=0,omega=1,f=0,eta=0,g=0):
        # get spatial dimension 
        ndim=self.getDim()
        # map kappa, omega, f, eta, g to the coefficients of the general PDE
        super(Helmholtz, self).setValue(A=kappa*kronecker(ndim),D=omega,Y=f,d=eta,y=g)
