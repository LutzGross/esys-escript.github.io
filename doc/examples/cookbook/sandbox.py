from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
from esys.pyvisi import Scene, DataCollector, Map, Camera 
from esys.pyvisi.constant import *
import os

#... set some parameters ...
xc=[0.0,0.0]
qc=50.e6
Tref=0.
rhocp=1
eta=10
kappa=1
tend=5.
# ... time, time step size and counter ...
t=0
h=0.1
i=0

#... generate domain ...
mydomain = Rectangle(l0=5,l1=0.00,n0=500, n1=1)
#... open PDE ...
mypde=LinearPDE(mydomain)
mypde.setSymmetryOn()
mypde.setValue(A=kappa*kronecker(mydomain),D=rhocp/h,d=eta,y=eta*Tref)
# ... set heat source: ....
x=mydomain.getX()
print x
#qH=qc*whereNegative(length(x-xc)-r)
