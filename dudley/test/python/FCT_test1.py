
##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2015 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

#
#  upwinding test moving a Gaussian hill around 
#
#     we solve U_,t + v_i u_,i =0 
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
from esys.dudley import Rectangle, Brick
from math import pi, ceil
NE=128
NE=4
DIM=2
THETA=0.5
OMEGA0=1.
ALPHA=pi/4
T0=0
T_END=2.*pi
dt=1e-3*10*10
E=1.e-3
TEST_SUPG=False or True


if DIM==2:
  dom=Rectangle(NE,NE)
else:
  dom=Brick(NE,NE,NE)
u0=dom.getX()[0]
# saveVTK("u.%s.vtu"%0,u=u0)
# print "XX"*80
dom.setX(2*dom.getX()-1)

# set initial value 
x=dom.getX()
r=sqrt(x[0]**2+(x[1]-1./3.)**2)
# u0=whereNegative(r-1./3.)*wherePositive(wherePositive(abs(x[0])-0.05)+wherePositive(x[1]-0.5))

x=Function(dom).getX()
if DIM == 2:
   V=OMEGA0*(x[0]*[0,-1]+x[1]*[1,0])
else:
   V=OMEGA0*(x[0]*[0,cos(ALPHA),0]+x[1]*[-cos(ALPHA),0,sin(ALPHA)]+x[2]*[0.,-sin(ALPHA),0.])
#===================
fc=TransportPDE(dom,num_equations=1,theta=THETA)
x=Function(dom).getX()
fc.setValue(M=Scalar(1.,Function(dom)),C=V)
#==============
if TEST_SUPG:
   supg=LinearSinglePDE(dom)
   supg.setValue(D=1.)
   supg.setSolverMethod(supg.LUMPING)
   dt_supg=inf(dom.getSize()/length(V))
   u_supg=u0*1.

c=0
# saveVTK("u.%s.vtu"%c,u=u0)
fc.setInitialSolution(u0)
t=T0
print("QUALITY FCT: time = %s pi"%(t/pi),inf(u0),sup(u0),integrate(u0))
while t<T_END:
    print("time step t=",t+dt)
    u=fc.solve(dt, verbose=True)
    print("QUALITY FCT: time = %s pi"%(t+dt/pi),inf(u),sup(u),integrate(u))
    if TEST_SUPG:
        #========== supg tests ================
        nn=max(ceil(dt/dt_supg),1.)
        dt2=dt/nn
        nnn=0
        while nnn<nn :
            supg.setValue(Y=u_supg+dt2/2*inner(V,grad(u_supg)))
            u2=supg.getSolution()
            supg.setValue(Y=u_supg+dt2*inner(V,grad(u2)))
            u_supg=supg.getSolution()
            nnn+=1
    c+=1
    t+=dt
    if TEST_SUPG: 
       print("QUALITY SUPG: time = %s pi"%(t/pi),inf(u_supg),sup(u_supg),integrate(u_supg))
       # saveVTK("u2.%s.vtu"%c,u=u,u_supg=u_supg)
    else:
       # saveVTK("u.%s.vtu"%c,u=u)
       pass
    # if c == 20: 1/0
