# $Id$
from mytools import *
from esys.escript import *
import esys.finley
#... set some parameters ...
x_c=[0.02,0.002]
r=0.001
q0=50.e6
Tref=0.
rhocp=2.6e6
eta=75.
kappa=240.
t_end=5.
# ...time step size and counter ...
h=0.1
i=0
t=0
#... generate domain ...
mydomain = esys.finley.Rectangle(l0=0.05,l1=0.01,n0=250, n1=50)
#... open PDE ...
mypde=Helmholtz(mydomain)
# ... set heat source: ....
x=mydomain.getX()
q=q0*(length(x-x_c)-r).whereNegative()
# ... set initial temperature ....
T=Tref
# ... start iteration:
while t<t_end:
      i+=1
      t+=h
      print "time step :",t
      mypde.setValue(kappa=kappa,omega=rhocp/h,f=q+rhocp/h*T,eta=eta,g=eta*Tref)
      T=mypde.getSolution()
      T.saveDX("T%d.dx"%i)
