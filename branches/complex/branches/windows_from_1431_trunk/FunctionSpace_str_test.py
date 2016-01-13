
import copy as cp

from esys.escript import *
from esys.finley import Rectangle
d=Rectangle(1,1)
y=Function(d)
#print type(y)
#x = cp.deepcopy(y)
#print type(x)
for i in range(50000):
   print (10*"="+" %d "+10*"=")%i
   sss=str(y)
   print sss

