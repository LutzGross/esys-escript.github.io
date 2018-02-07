
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

#
#  upwinding test moving a Gaussian hill around 
#
#     we solve U_,t - E *u_,ii + v_i u_,i =0 (E is small)
#
#  the solution is given as   u(x,t)=1/(4*pi*E*t)^{dim/2} * exp ( - |x-x_0(t)|^2/(4*E*t) ) 
#
#   where x_0(t) = [ cos(OMEGA0*T0)*0.5,-sin(OMEGA0*T0)*0.5 ] and v=[-y,x]*OMEGA0 for dim=2 and
#
#         x_0(t) = [ cos(OMEGA0*T0)*0.5,-sin(OMEGA0*T0)*0.5 ] and v=[-y,x]*OMEGA0 for dim=3
#
#  the solution is started from some time T0>0.
#
#  We are using five quality messurements for u_h
#
#     - inf(u_h) > 0
#     - sup(u_h)/sup(u(x,t)) = sup(u_h)*(4*pi*E*t)^{dim/2} ~ 1 
#     - integrate(u_h) ~ 1
#     - | x_0h-x_0 | ~ 0    where x_0h = integrate(x*u_h)
#     - sigma_h/4*E*t ~ 1 where sigma_h=sqrt(integrate(length(x-x0h)**2 * u_h) * (DIM==3 ? sqrt(2./3.) :1 )
#
#
from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, TransportPDE
from esys.finley import Rectangle, Brick
from esys.weipa import saveVTK
from math import pi, ceil
NE=128
NE=64
DIM=2
THETA=0.5
OMEGA0=1.
ALPHA=pi/4
T0=0.5*pi
T_END=2.5*pi
dt=1e-3*10
E=1.e-3
TEST_SUPG=False or True


def getCenter(t):
   if DIM==2:
      x0=[cos(OMEGA0*t)*0.5,-sin(OMEGA0*t)*0.5]
      x0=[-sin(OMEGA0*t)*0.5,cos(OMEGA0*t)*0.5]
   else:
      x0=[cos(ALPHA)*cos(OMEGA0*t)*0.5,-sin(OMEGA0*t)*0.5,-sin(ALPHA)*cos(OMEGA0*t)*0.5]
   return x0
def QUALITY(t,u_h):
     dom=u_h.getDomain()
     x=dom.getX()
     a=inf(u_h)
     b=sup(u_h)*(4*pi*E*t)**(DIM/2.)-1.
     c=integrate(u_h,Function(dom))-1.
     x0=getCenter(t)
     x0h=integrate(x*u_h,Function(dom))
     d=length(x0-x0h)
     sigma_h2=sqrt(integrate(length(x-x0h)**2 * u_h, Function(dom)))
     if DIM == 3: sigma_h2*=sqrt(2./3.)
     e=sigma_h2/sqrt(4*E*t)-1             
     # return a,b,c,e,1./(4*pi*E*t)**(DIM/2.)
     return d,e
     # return a,b,c,d,e
      



if DIM==2:
  dom=Rectangle(NE,NE)
else:
  dom=Brick(NE,NE,NE)
dom.setX(2*dom.getX()-1)

# set initial value 
x=dom.getX()
u0=1/(4.*pi*E*T0)**(DIM/2.)*exp(-length(dom.getX()-getCenter(T0))**2/(4.*E*T0)) 

print("QUALITY ",QUALITY(T0,u0))

x=Function(dom).getX()
if DIM == 2:
   V=OMEGA0*(x[0]*[0,-1]+x[1]*[1,0])
else:
   V=OMEGA0*(x[0]*[0,cos(ALPHA),0]+x[1]*[-cos(ALPHA),0,sin(ALPHA)]+x[2]*[0.,-sin(ALPHA),0.])
#===================
fc=TransportPDE(dom,num_equations=1,theta=THETA)
x=Function(dom).getX()
fc.setValue(M=Scalar(1.,Function(dom)),C=V,A=-Scalar(E,Function(dom))*kronecker(dom))
#==============
if TEST_SUPG:
   supg=LinearSinglePDE(dom)
   supg.setValue(D=1.)
   supg.setSolverMethod(supg.LUMPING)
   dt_supg=1./(1./inf(dom.getSize()/length(V))+1./inf(dom.getSize()**2/E))*0.3
   u_supg=u0*1.

c=0
saveVTK("u.%s.vtu"%c,u=u0)
fc.setInitialSolution(u0)
t=T0
while t<T_END:
    print("time step t=",t+dt)
    u=fc.solve(dt)
    if TEST_SUPG:
        #========== supg tests ================
        nn=max(ceil(dt/dt_supg),1.)
        dt2=dt/nn
        nnn=0
        while nnn<nn :
            supg.setValue(X=-dt2/2*E*grad(u_supg),Y=u_supg+dt2/2*inner(V,grad(u_supg)))
            u2=supg.getSolution()
            supg.setValue(X=-dt2*E*grad(u2),Y=u_supg+dt2*inner(V,grad(u2)))
            u_supg=supg.getSolution()
            nnn+=1
    c+=1
    t+=dt
    print("QUALITY FCT: time = %s pi"%(t/pi),QUALITY(t,u), end=' ')
    if TEST_SUPG: 
       print("QUALITY SUPG: ",QUALITY(t,u_supg))
       # saveVTK("u.%s.vtu"%c,u=u,u_supg=u_supg)
    else:
       # saveVTK("u.%s.vtu"%c,u=u)
       pass
    # if c == 20: 1/0
