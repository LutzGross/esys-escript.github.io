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
"""
Damage mechanics 
"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDESystem
from esys.finley import Rectangle, Brick
from esys.weipa import saveVTK
from esys.escript import unitsSI as U
from math import pi, ceil
import sys
import time

# ======================= Default Values ==================================================
H=1.*U.m                     # height
L=1.*H                           # length
NE_H=1                          # number of elements in H-direction. 
NE_L=int(ceil(L*NE_H/H))
CASE=1

#Boundary conditions: 
#   axial loading: they applied a stress inversely proportional to the acoustic emission rate. We could have the axial forcing a stress or velocity inversely proportional to dalpha/dt (only when it is positive, and with the applied forcing rate going to zero when damage accumulation rate goes to a value we can determine in a test run with constant forcing). If this is to challenging or time consuming we could have a constant axial strain rate with very short time steps (at least when alpha increases above 0.3).

#Variables calculated and written to an output file:
#   time
#   differential stress (S_33-S_11)
#   deviatoric stress (S_33 - p)
#   Axial and transverse strain
#   damage and damage rate


T_END=60000000.0*U.sec                       # end time

if CASE==1 or CASE==2:
   C_V = 2.e-5/(U.Mega*U.Pa)  *0
   C_D = 3/U.sec   
   C_1 = 1e-12/U.sec   
   C_2 = 0.03
   XI_0 = -0.56
   LAME_0 = 29e9*U.Pa
   MU_0 = 19e9*U.Pa
   if CASE==2:
      ALPHA_0 = 0.999
   else:
      ALPHA_0 = 0.0
   RHO = 2800*U.kg/U.m**3
   G = 10*U.m/U.sec**2     *0
   SIGMA_N=-50.*U.Mega*U.Pa
   DIM=3                          
   VMAX=-1.*U.m/U.sec/166667.
   VMAX=-1.*U.m/U.sec/100000.
   DT_MAX=50.*U.sec
   DT=DT_MAX/1000000.
   xc=[L/2,L/2,H/2]
   WWW=min(H,L)*0.01

else:
   C_V = 3e-11/U.Pa 
   C_D = 5/U.sec
   C_1 = 1e-12/U.sec 
   C_2 = 0.03
   XI_0 = -0.8
   LAME_0 = 46e9*U.Pa
   MU_0 = 30e9*U.Pa
   if CASE==3:
      ALPHA_0 = 0.01
   else:
      ALPHA_0 = 0.0
   RHO = 2800*U.kg/U.m**3
   G = 10*U.m/U.sec**2     *0
   VMAX=-1*U.m/U.sec
   DT_MAX=500.*U.sec
   DT=DT_MAX/100000.
   SIGMA_N=0
   DIM=2                      
   xc=[L/2,H/2]
   WWW=min(H,L)*0.08

VERBOSE=True
DT_VIS=T_END/100                # time distane between two visulaization files
DN_VIS=1                        # maximum counter increment between two visulaization files
VIS_DIR="results"               # name of the director for vis files

ODE_TOL=0.01
ODE_ITER_TOL=1.e-8
ODE_ITER_MAX=15
DEPS_MAX=0.01
TOL_DU=1e-8
UPDATE_OPERATOR = False


diagnose=FileWriter("diagnose.csv",append=False)
#===================================
S=0.5*XI_0*((2.*MU_0+3.*LAME_0)/(3.-XI_0**2) + LAME_0)
GAMMA_M=S + sqrt(S**2+2.*MU_0*(2.*MU_0+3.*LAME_0)/(3.-XI_0**2))
def solveODE(u0, a, b, dt):
    """
    solves du/dt=a*exp(b*u)   u(t=0)=u0 and return approximation at t=dt

    we sue backwards Euler by solving     u-u0 = dt * a*exp(b*u)
    with newton scheme                    u <- u - (u-u0 - dt * a*exp(b*u)) / (1-dt*a*b*exp(b*u))
    """
    u=u0.copy()
    norm_du=1.
    norm_u=0.
    n=0
    while norm_du > ODE_ITER_TOL * norm_u:
         H=-dt*a*exp(b*u)
         du=-(u-u0+H)/(1+b*H)
         u+=du
         norm_du = Lsup(du)
         norm_u = Lsup(u)
         n+=1
         if n>ODE_ITER_MAX: raise ValueError("ODE iteration failed.")
    print("\tODE iteration completed after %d steps."%(n,))
    return u

#======================
t=0         # time stamp
n=0         # time step counter
dt=DT       # current time step size
t_vis=0     
n_vis=0
counter_vis=0
mkDir(VIS_DIR)
#=========================
#
#   set up domain 
#
if DIM==2:
    dom=Rectangle(NE_L,NE_H,l0=L,l1=H,order=1,optimize=True)
else:
    dom=Brick(NE_L,NE_L,NE_H,l0=L,l1=L,l2=H,order=1,optimize=True)

BBOX=boundingBox(dom)
DIM=dom.getDim()
x=dom.getX()
#
#   initial values:
#
sigma=Tensor(0.,Function(dom))
eps_e=Tensor(0.,Function(dom))
if CASE==2 or CASE==3:
  alpha=ALPHA_0*exp(-length(Function(dom).getX()-xc)**2/WWW**2)
else:
  alpha=Scalar(ALPHA_0,Function(dom))

pde=LinearPDESystem(dom)
pde.setSymmetryOn()
pde.getSolverOptions().setSolverMethod(pde.getSolverOptions().DIRECT)



fixed_v_mask=Vector(0,Solution(dom))
v0=Vector(0.,ContinuousFunction(dom))
if CASE == 1 or CASE==2:
    for d in range(DIM):
       fixed_v_mask+=whereZero(x[d]-BBOX[d][0])*unitVector(d,DIM)
       if d == DIM-1:
          fixed_v_mask+=whereZero(x[d]-BBOX[d][1])*unitVector(d,DIM)
          v0[d]=(x[d]-BBOX[d][0])/(BBOX[d][1]-BBOX[d][0])*VMAX
else:
    for d in range(DIM):
        fixed_v_mask+=whereZero(x[d]-BBOX[d][0])*unitVector(d,DIM)
        if d == 0:
           fixed_v_mask+=whereZero(x[d]-BBOX[d][1])*unitVector(d,DIM)
           v0[d]=(x[d]-BBOX[d][0])/(BBOX[d][1]-BBOX[d][0])*VMAX

pde.setValue(Y=-G*RHO*kronecker(DIM)[DIM-1], q=fixed_v_mask)
du=Vector(0.,Solution(dom))
u=Vector(0.,Solution(dom))
norm_du=0.
deps=Tensor(0,Function(dom))
i_eta=0
#
#  let the show begin:
#
k3=kronecker(DIM)
k3Xk3=outer(k3,k3)
alpha_old=alpha
dt_old=None
diagnose.write("t, -e22, e11, s00-s22, mu_eff, lame_eff, xi, gamma, alpha, alpha_dot\n")

while t<T_END:

    print("start time step %d"%(n+1,))

    eps_e_old = eps_e
    sigma_old = sigma 
    alpha_old, alpha_oold = alpha, alpha_old
    #  start the iteration for deps on a time step: deps from the last time step is used as an initial guess:
    iter=0
    norm_ddu=norm_du
    while norm_ddu > TOL_DU * norm_du or iter == 0 :

        print("\t start iteration step %d:"%iter)
        eps_e = eps_e_old + deps-(dt/2)*i_eta*deviatoric(sigma)
       
        I1=trace(eps_e)
        sqI2=length(eps_e)
        xi=safeDiv(I1,sqI2)
        i_xi=safeDiv(sqI2,I1)
        # update damage parameter:
        m=wherePositive(xi-XI_0)
        a=sqI2**2*(xi-XI_0)*(m*C_D + (1-m)* C_1)
        b=(1-m)*(1./C_2)

        alpha=solveODE(alpha_old, a,b, dt)
        alpha_dot=(alpha-alpha_old)/dt
        i_eta = clip(2*C_V*alpha_dot,minval=0.)

        if inf(alpha) < -EPSILON*10:
            raise ValueError("alpha<0")
        if sup(alpha)  > 1:
            raise ValueError("alpha > 1")
        # step size for the next time step:

        gamma=alpha*GAMMA_M
        lame=LAME_0
        mu=MU_0*(1-alpha)

        lame_eff=lame-gamma*i_xi
        mu_eff=mu-gamma*xi/2.

        print("\talpha = [ %e, %e]"%(inf(alpha),sup(alpha)))
        print("\tmu_eff = [ %e, %e]"%(inf(mu_eff),sup(mu_eff)))
        print("\tlame_eff = [ %e, %e]"%(inf(lame_eff),sup(lame_eff)))
        print("\txi = [ %e, %e]"%(inf(xi),sup(xi)))
        print("\tgamma = [ %e, %e]"%(inf(gamma),sup(gamma)))

        if inf(mu_eff) < 0:
            raise ValueError("mu_eff<0")

        sigma = 2*mu_eff*eps_e+lame_eff*trace(eps_e)*k3 

        if (UPDATE_OPERATOR) :
            pde.setValue(A = mu_eff * ( swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3) ) + lame_eff*k3Xk3)
        else:
            pde.setValue(A = mu * ( swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3) ) + lame*k3Xk3)

        pde.setValue(X=-sigma, y=SIGMA_N*dom.getNormal(), r=dt*v0-du)
        

        ddu=pde.getSolution()
        deps+=symmetric(grad(ddu))
        du+=ddu
        norm_ddu=Lsup(ddu)
        norm_du=Lsup(du)
        print("\t displacement change update = %e of %e"%(norm_ddu, norm_du))
        iter+=1
    print("deps =", inf(deps),sup(deps))
    u+=du
    n+=1
    t+=dt
    #=========== this is a test for triaxial test ===========================
    print("\tYYY t = ", t)
    a =(SIGMA_N-lame_eff*VMAX*t)/(lame_eff+mu_eff)/2
    #=========== this is a test for triaxial test ===========================
    print("\tYYY a = ", meanValue(a))
    print("\tYYY eps00 = ",meanValue( eps_e[0,0]))
    print("\tYYY eps11 = ",meanValue( eps_e[1,1]))
    print("\tYYY eps22 = num/exact", meanValue(eps_e[2,2]), VMAX*t)
    print("\tYYY eps_kk = num/exact", meanValue(trace(eps_e)), meanValue(VMAX*t+2*a))
    print("\tYYY sigma11 = num/exact", meanValue(sigma[1,1]), meanValue(lame_eff*(VMAX*t+2*a)+2*mu_eff*a))
    print("\tYYY sigma22 = num/exact", meanValue(sigma[2,2]), meanValue(lame_eff*(VMAX*t+2*a)+2*mu_eff*VMAX*t))
    print("\tYYY linear Elastic equivalent num/exact=",meanValue(sigma[2,2]-sigma[0,0]-(sigma_old[2,2]-sigma_old[0,0]))/meanValue(eps_e[2,2]-eps_e_old[2,2]), meanValue(mu_eff*(3*lame_eff+2*mu_eff)/(lame_eff+mu_eff)))
    diagnose.write(("%e,"*10+"\n")%(t, meanValue(-eps_e[2,2]), 
                                      meanValue(eps_e[1,1]),
                                      meanValue(sigma[0,0]-sigma[2,2]), 
                                      meanValue(mu_eff),
                                      meanValue(lame_eff),
                                      meanValue(xi),
                                      meanValue(gamma),
                                      meanValue(alpha),
                                      meanValue(alpha_dot)))
    print("time step %s (t=%s) completed."%(n,t))
    #
    #  .... visualization
    #
    if t>=t_vis or n>n_vis:
      saveVTK(os.path.join(VIS_DIR,"state.%d.vtu"%counter_vis),u=u, dalpha=alpha, I1=trace(eps_e), I2=length(eps_e)**2, xi=safeDiv(trace(eps_e),length(eps_e)))
      print("visualization file %d for time step %e generated."%(counter_vis,t))
      counter_vis+=1
      t_vis+=DT_VIS
      n_vis+=DN_VIS
    # 
    #   control time step size:
    #
    ss=sup(length(deps))
    if ss>0:
       dt_new=DEPS_MAX/ss*dt
       print("\ttime step size to control strain increment %s."%(dt_new,))
    else:
       dt_new=dt
    if dt_old != None:
          dd_alpha=2.*dt_old*(alpha-alpha_old)+(alpha_oold-alpha_old)*dt/(dt*dt_old*(dt_old+dt))
          norm_alpha=Lsup(alpha)
          fac=Lsup(dd_alpha)
          if norm_alpha > 0: fac*=1./norm_alpha
          if fac>0:
             error=fac*0.5*dt**2
             print("\testimated local error for time step size %e is %e"%(dt,error))
             dt_new=min(dt_new,sqrt(ODE_TOL*2/fac) )
          else:
             dt_new=DT_MAX

    dt_new=min(max(dt_new,dt/5),dt*5,DT_MAX) # aviod rapit changes
    print("\tINFO: new time step size %e"%dt_new)
    dt, dt_old=dt_new, dt
    # dom.setX(dom.getX()+du)
