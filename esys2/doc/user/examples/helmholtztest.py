# $Id$
from mytools import Helmholtz
from esys.escript import Lsup
from esys.finley import Rectangle
#... set some parameters ...
kappa=1.
omega=0.1
eta=10.
#... generate domain ...
mydomain = Rectangle(l0=5.,l1=1.,n0=50, n1=10)
#... open PDE and set coefficients ...
mypde=Helmholtz(mydomain)
n=mydomain.getNormal()
x=mydomain.getX()
mypde.setValue(kappa,omega,omega*x[0],eta,kappa*n[0]+eta*x[0])
#... calculate error of the PDE solution ...
u=mypde.getSolution()
print "error is ",Lsup(u-x[0])
# output should be similar to "error is 1.e-7" 
