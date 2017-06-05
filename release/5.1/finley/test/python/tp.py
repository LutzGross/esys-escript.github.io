
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, TransportPDE
from esys.finley import Rectangle
from esys.weipa import saveVTK

# dom=Rectangle(12,8,l0=1.5)
# dom=Rectangle(24,16,l0=1.5)
dom=Rectangle(48,32,l0=1.5)
saveDataCSV("t.csv",x=dom.getX(), rho=length(dom.getX()))
1/0
# dom=Rectangle(8*48,8*32,l0=1.5)
# dom=Rectangle(120,80,l0=1.5)
V=Scalar(1.,Function(dom))*[-1.,0]
THETA=0.
fc=TransportPDE(dom,num_equations=1,theta=THETA)
fc.setTolerance(1.e-12)
fc.setValue(M=Scalar(1.,Function(dom)),C=V)
x=dom.getX()
x_0=[0.5,0.5]
sigma=0.075
u0=1.
for i in range(dom.getDim()):
    u0=u0*exp(-(x[i]-x_0[i])**2/sigma**2)

u0=whereNonPositive(abs(x[0]-0.4)-0.2)*whereNonPositive(abs(x[1]-0.5)-0.2)
# f1=0.5
# f2=2.
# u0=f2*clip(x[0]-0.5,0.)-clip(0.5-x[0],0.)*f1+f1*0.5
# u0=exp(-3*(x[0]-2.)**2)
# u0=x[0]
u0/=Lsup(u0)
c=0
saveVTK("u.%s.vtu"%c,u=u0)
fc.setInitialSolution(u0)

t_end=0.6
dt=2.49999e-2*0+6.2499999e-02/4
dt_out=2.49999e-2*0+6.2499999e-02/4
c_stop=1
n_out=int(t_end/dt+0.5)
print(n_out)
t=0.
t_out=0
c_out=0
c=0
print(t,": range u",inf(u0),sup(u0),integrate(u0,Function(dom)))
while t<t_end and c< c_stop:
    print("time step t=",t+dt)
    u=fc.solve(dt)
    print(t+dt,": range u",inf(u),sup(u),integrate(u,Function(dom)))
    c+=1
    t+=dt
    if t>=t_out+dt_out:
         c_out,t_out=c_out+1,t_out+dt_out
         saveVTK("u.%s.vtu"%c_out,u=u)
         print("write time step ",c,"(t=%s) to file u.%s.vtu"%(t,c_out))

if True:
   pde=LinearPDE(dom)
   pde.setValue(D=1.,C=-THETA*dt*V)
   pde.setTolerance(1e-12)
   t=0.
   t_out=0
   c_out=0
   c=0
   u=u0
   print(t,": range u2",inf(u0),sup(u0),integrate(u0,Function(dom)))
   while t<t_end and c< c_stop:
       print("time step t=",t+dt)
       pde.setValue(Y=u+(1.-THETA)*dt*inner(V,grad(u)))
       u=pde.getSolution(verbose=True)
       print(t+dt,": range u2",inf(u),sup(u),integrate(u,Function(dom)))
       c+=1
       t+=dt
       if t>=t_out+dt_out:
         c_out,t_out=c_out+1,t_out+dt_out
         saveVTK("u2.%s.vtu"%c_out,u=u)
         print("write time step ",c,"(t=%s) to file u2.%s.vtu"%(t,c_out))
