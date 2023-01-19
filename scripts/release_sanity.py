
# Very basic sanity checks following a build
# It is important that this script does not create any files unless
# it is _certain_ they are removed when finished or failed.
# We do not want the source directory polluted by actions here


from esys.escript import *
from esys.escript.linearPDEs import Poisson
import esys.ripley as ripley
import esys.finley as finley
import esys.oxley as oxley
import esys.speckley as speckley
from esys.weipa import saveVTK

mydomain = finley.Rectangle(l0=1.,l1=1.,n0=40, n1=20)
# define characteristic function of Gamma^D
x = mydomain.getX()
gammaD = whereZero(x[0])+whereZero(x[1])
# define PDE and get its solution u
mypde = Poisson(domain=mydomain)
mypde.setValue(f=1,q=gammaD)
u = mypde.getSolution()
# write u to an external file  
#saveVTK("u.vtu",sol=u)

print("Passed the release_sanity test!")


