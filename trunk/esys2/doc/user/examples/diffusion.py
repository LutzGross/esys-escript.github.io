# $Id$
from mytools import Helmholtz
from esys.escript import Lsup
from esys.finley import Rectangle
#... set some parameters ...
xc=[0.02,0.002]
r=0.001
qc=50.e6
Tref=0.
rhocp=2.6e6
eta=75.
kappa=240.
tend=5.
# ... time, time step size and counter ...
t=0
h=0.1
i=0
#... generate domain ...
mydomain = Rectangle(l0=0.05,l1=0.01,n0=250, n1=50)
#... open PDE ...
mypde=Helmholtz(mydomain)
# ... set heat source: ....
x=mydomain.getX()
q=qc*(length(x-xc)-r).whereNegative()
# ... set initial temperature ....
T=Tref
# ... start iteration:
while t<tend:
      i+=1
      t+=h
      print "time step :",t
      mypde.setValue(kappa=kappa,omega=rhocp/h,f=q+rhocp/h*T,eta=eta,g=eta*Tref)
      T=mypde.getSolution()
      T.saveDX("T%d.dx"%i)