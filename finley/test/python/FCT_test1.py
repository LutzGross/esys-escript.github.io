
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
from esys.escript.linearPDEs import TransportPDE, SolverOptions
from esys.finley import Rectangle, Brick
#from esys.ripley import Rectangle, Brick
from esys.weipa import saveVTK
from math import pi, ceil
NE=128
#NE=4
DIM=2
THETA=0.5
OMEGA0=1.
ALPHA=pi/4
T0=0
T_END=2.*pi
dt=1e-3*10*10
E=1.e-3


dom=Rectangle(NE,NE)
u0=dom.getX()[0]
# saveVTK("u.%s.vtu"%0,u=u0)
# print "XX"*80

# set initial value 
#dom.setX(2*dom.getX()-1)
#x=dom.getX()
#r=sqrt(x[0]**2+(x[1]-1./3.)**2)
#u0=whereNegative(r-1./3.)*wherePositive(wherePositive(abs(x[0])-0.05)+wherePositive(x[1]-0.5))

#x=Function(dom).getX()
#if DIM == 2:
#   V=OMEGA0*(x[0]*[0,-1]+x[1]*[1,0])
#else:
#   V=OMEGA0*(x[0]*[0,cos(ALPHA),0]+x[1]*[-cos(ALPHA),0,sin(ALPHA)]+x[2]*[0.,-sin(ALPHA),0.])

x=dom.getX()

R0=0.15
#cylinder:
X0=0.5
Y0=0.75
r=sqrt((x[0]-X0)**2+(x[1]-Y0)**2)/R0
u0=whereNegative(r-1)*wherePositive(wherePositive(abs(x[0]-X0)-0.025)+wherePositive(x[1]-0.85))
# cone:
X0=0.5
Y0=0.25
r=sqrt((x[0]-X0)**2+(x[1]-Y0)**2)/R0
u0=u0+wherePositive(1-r)*(1-r)
#hump
X0=0.25
Y0=0.5
r=sqrt((x[0]-X0)**2+(x[1]-Y0)**2)/R0
u0=u0+1./4.*(1+cos(pi*clip(r,maxval=1)))

x=Function(dom).getX()
V=OMEGA0*((0.5-x[0])*[0,1]+(0.5-x[1])*[-1,0])
#===================

fc=TransportPDE(dom,numEquations=1)
fc.getSolverOptions().setVerbosityOn()
#fc.getSolverOptions().setODESolver(SolverOptions.BACKWARD_EULER)
fc.getSolverOptions().setODESolver(SolverOptions.LINEAR_CRANK_NICOLSON)
fc.getSolverOptions().setODESolver(SolverOptions.CRANK_NICOLSON)
x=Function(dom).getX()
fc.setValue(M=1,C=V)

c=0
saveVTK("u.%s.vtu"%c,u=u0)
fc.setInitialSolution(u0)
dt=fc.getSafeTimeStepSize() 
#dt=1.e-3
print("dt = ",dt)
t=T0
print("QUALITY FCT: time = %s pi"%(t/pi),inf(u0),sup(u0),integrate(u0))
#T_END=200*dt
while t<T_END:
   
    print("time step t=",t+dt)
    u=fc.getSolution(dt)
    print("QUALITY FCT: time = %s pi"%(t+dt/pi),inf(u),sup(u),integrate(u))
    saveVTK("u.%s.vtu"%(c+1,),u=u)
    c+=1
    t+=dt
