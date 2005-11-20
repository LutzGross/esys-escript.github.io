# $Id$
from esys.linearPDEs import LinearPDE
import numarray
class Helmholtz(LinearPDE):
   def setValue(self,kappa=0,omega=1,f=0,eta=0,g=0):
        # get spatial dimension 
        ndim=self.getDim()
        # get kronecker symbol from numarray: kronecker[i,j]=1 for i=j and =0 else
        kronecker=numarray.identity(ndim)
        # map kappa, omega, f, eta, g to the coefficients of the general PDE
        self._LinearPDE__setValue(A=kappa*kronecker,D=omega,Y=f,d=eta,y=g)
