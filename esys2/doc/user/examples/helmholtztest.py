# $Id$
from mytools import *
from esys.escript import *
import esys.finley
#... set some parameters ...
omega=0.1
eta=10.
#... generate domain ...
mydomain = esys.finley.Rectangle(l0=5.,l1=1.,n0=50, n1=10)
#... open PDE and set coefficients ...
mypde=Helmholtz(mydomain)
n=mydomain.getNormal()
x=mydomain.getX()
mypde.setValue(1,omega,omega*x[0],eta,n[0]+eta*x[0])
#... calculate error of the PDE solution ...
u=mypde.getSolution()
print "error is ",Lsup(u-x[0])
