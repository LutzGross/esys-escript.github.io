# $Id$
import esys.finley
from esys.linearPDEs import Poisson
# generate domain:
mydomain = esys.finley.Rectangle(l0=1.,l1=1.,n0=40, n1=20)
# define characteristic function of Gamma^D
x = mydomain.getX()
gammaD = x[0].whereZero()+x[1].whereZero()
# define PDE and get its solution u
mypde = Poisson(domain=mydomain,f=1,q=gammaD)
u = mypde.getSolution()
# write u to an external file
u.saveDX("u.dx")
