########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################
"""
Damage mechanics 
"""
__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDESystem
from esys.finley import Rectangle, Brick
from esys.escript import unitsSI as U
from math import pi, ceil
import sys
import time

# ======================= Default Values ==================================================
H=1.*U.m                     # height
L=1.*H                           # length
NE_H=5                           # number of elements in H-direction. 
NE_L=int(ceil(L*NE_H/H))

#Boundary conditions: 
#	axial loading: they applied a stress inversely proportional to the acoustic emission rate. We could have the axial forcing a stress or velocity inversely proportional to dalpha/dt (only when it is positive, and with the applied forcing rate going to zero when damage accumulation rate goes to a value we can determine in a test run with constant forcing). If this is to challenging or time consuming we could have a constant axial strain rate with very short time steps (at least when alpha increases above 0.3).

#Variables calculated and written to an output file:
#	time
#	differential stress (S_33-S_11)
#	deviatoric stress (S_33 - p)
#	Axial and transverse strain
#	damage and damage rate



if True:
   C_V = 1.e-5/(U.Mega*U.Pa)
   C_D = 3/U.sec
   C_1 = 1e-12/U.sec
   C_2 = 0.03
   XI_0 = -0.56
   LAME_0 = 29e9*U.Pa
   MU_0 = 19e9*U.Pa
   ALPHA_0 = 0.0
   RHO = 2800*U.kg/U.m**3
   G = 10*U.m/U.sec**2     *0
   DT=1.*U.sec
   SIGMA_N=50.*U.Mega*U.Pa
   DIM=3                          
   CASE=1
   VMAX=-0.1*U.m/U.sec
   DT_MAX=10.*U.sec

else:
   C_V = 3e-11/U.Pa 
   C_D = 5/U.sec
   C_1 = 1e-12/U.sec
   C_2 = 0.03
   XI_0 = -0.8
   LAME_0 = 46e9*U.Pa
   MU_0 = 30e9*U.Pa
   ALPHA_0 = 0.1
   RHO = 2800*U.kg/U.m**3
   G = 10*U.m/U.sec**2     *0
   DT=1.*U.sec/100.
   VMAX=-0.01*U.m/U.sec
   DT_MAX=10.*U.sec
   SIGMA_N=0
   DIM=2                      
   CASE=0

T_END=15.*U.sec                       # end time
VERBOSE=True
DT_VIS=T_END/100                # time distane between two visulaization files
DN_VIS=5                        # maximum counter increment between two visulaization files
VIS_DIR="results"               # name of the director for vis files

ODE_TOL=0.01
ODE_ITER_TOL=1.e-8
ODE_ITER_MAX=15

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
         if n>ODE_ITER_MAX: raise ValueError,"ODE iteration failed."
    print "\tODE iteration completed after %d steps."%(n,)
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
    dom=Rectangle(NE_L,NE_H,l0=L,l1=H,order=-1,optimize=True)
else:
    dom=Brick(NE_L,NE_L,NE_H,l0=L,l1=L,l2=H,order=-1,optimize=True)

BBOX=boundingBox(dom)
DIM=dom.getDim()
x=dom.getX()
#
#   initial values:
#
sigma=Tensor(0.,Function(dom))
eps_e=Tensor(0.,Function(dom))
alpha=Scalar(ALPHA_0,Function(dom))
u=Vector(0.,ContinuousFunction(dom))

pde=LinearPDESystem(dom)
pde.setSymmetryOn()
pde.getSolverOptions().setSolverMethod(pde.getSolverOptions().DIRECT)



fixed_v_mask=Vector(0,Solution(dom))
u0=Vector(0.,ContinuousFunction(dom))
if CASE == 1:
    for d in range(DIM):
       fixed_v_mask+=whereZero(x[d]-BBOX[d][0])*unitVector(d,DIM)
       if d == DIM-1:
          fixed_v_mask+=whereZero(x[d]-BBOX[d][1])*unitVector(d,DIM)
          u0[d]=(x[d]-BBOX[d][0])/(BBOX[d][1]-BBOX[d][0])*VMAX
else:
    for d in range(DIM):
        fixed_v_mask+=whereZero(x[d]-BBOX[d][0])*unitVector(d,DIM)
        if d == 0:
           fixed_v_mask+=whereZero(x[d]-BBOX[d][1])*unitVector(d,DIM)
           u0[d]=(x[d]-BBOX[d][0])/(BBOX[d][1]-BBOX[d][0])*VMAX

pde.setValue(Y=-G*RHO*kronecker(DIM)[DIM-1], q=fixed_v_mask, r=u0, y=SIGMA_N*dom.getNormal())
#
#  let the show begin:
#
k3=kronecker(DIM)
k3Xk3=outer(k3,k3)
alpha_old=alpha
dt_old=None
if CASE == 1:
   diagnose.write("t,s22-s00, s22-p, e00, e11, alpha, alpha_dot\n")
else:
   diagnose.write("t, eps00, eps11, error0, sigma11/tau\n")
while t<T_END:

    print "start time step %d"%(n+1,)

    I1=trace(eps_e)
    sqI2=length(eps_e)
    print "II2 =",sqI2
    xi=safeDiv(I1,sqI2)
    i_xi=safeDiv(sqI2,I1)
    # update damage parameter:
    m=wherePositive(xi-XI_0)
    a=sqI2**2*(xi-XI_0)*(m*C_D + (1-m)* C_1)
    b=(1-m)*(1./C_2)

    alpha, alpha_old, alpha_oold =solveODE(alpha, a,b, dt), alpha, alpha_old
    alpha_dot=(alpha-alpha_old)/dt


    # step size for the next time step:

    print "\tXXX alpha = [ %e, %e]"%(inf(alpha),sup(alpha))
    if sup(alpha)  >1:
        raise ValueError,"damage parameter > 1 detected."

    gamma=alpha*GAMMA_M
    lame=LAME_0
    mu=MU_0*(1-alpha)

    lame_eff=lame-gamma*i_xi
    mu_eff=mu-gamma*xi
    print "\tmu_eff = [ %e, %e]"%(inf(mu_eff),sup(mu_eff))
    print "\tlame_eff = [ %e, %e]"%(inf(lame_eff),sup(lame_eff))

    if CASE == 1:
      diagnose.write(("%e,"*7+"\n")%(t,meanValue(sigma[2,2]-sigma[0,0]), meanValue(deviatoric(sigma)[2,2]), meanValue(eps_e[0,0]), meanValue(eps_e[1,1]),
                    meanValue(alpha), meanValue(alpha_dot)))
    else:
       if n>1:
          print "\t xi = [ %e, %e]"%(inf(xi),sup(xi))
          print "\t eps00 = ", eps_e[0,0]
          print "\t eps10 = ", eps_e[1,0]
          print "\t eps01 = ", eps_e[0,1]
          print "\t eps11 = ", eps_e[1,1]
          print "\t eps11/eps00 = ", eps_e[1,1]/eps_e[0,0]
          print" error:",Lsup(((2*mu-gamma*safeDiv(I1,sqI2))*eps_e[1,1]-(gamma*sqI2-lame*I1)) )/Lsup(length(sigma))
          print" error:",Lsup(sigma[1,1])/Lsup(length(sigma))
          diagnose.write("%e, %e, %e, %e, %e\n"%(t,inf(eps_e[0,0]),inf(eps_e[1,1]), Lsup(((2*mu-gamma*safeDiv(I1,sqI2))*eps_e[1,1]-(gamma*sqI2-lame*I1)) )/Lsup(length(sigma)), Lsup(sigma[1,1])/Lsup(length(sigma))))

    if inf(mu_eff) <= 0: 
        raise ValueError,"Material failed due to non-positive mu_eff."
    if inf(lame_eff) <= 0: 
        raise ValueError,"Material failed due to non-positive lame_eff."

    i_eta = clip(2*C_V*alpha_dot,minval=0.)

    # update strain:

    H = safeDiv(I1,sqI2**2)*eps_e-k3
    pde.setValue(A = mu_eff * ( swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3) ) + lame_eff*k3Xk3 + gamma*i_xi * outer(H,H),
                 X=- (sigma + dt * mu_eff * i_eta * gamma * safeDiv(length(deviatoric(eps_e))**2,sqI2) * H ) ) 

    du=pde.getSolution()*dt

    eps_e=eps_e+symmetric(grad(du))-dt/2*i_eta*deviatoric(sigma)
    sigma=2*mu_eff*eps_e+lame_eff*trace(eps_e)*k3
    print "time step %s (t=%s) completed."%(n,t)
    n+=1
    t+=dt
    #
    #  .... visualization
    #
    if t>=t_vis or n>n_vis:
      saveVTK(os.path.join(VIS_DIR,"state.%d.vtu"%counter_vis),alpha=alpha, I1=trace(eps_e), I2=length(eps_e)**2)
      print "visualization file %d for time step %e generated."%(counter_vis,t)
      counter_vis+=1
      t_vis+=DT_VIS
      n_vis+=DN_VIS
    # 
    #   control time step size:
    #
    if dt_old == None:
          dt_new=dt
    else:
          dd_alpha=2.*dt_old*(alpha-alpha_old)+(alpha_oold-alpha_old)*dt/(dt*dt_old*(dt_old+dt))
          norm_alpha=Lsup(alpha)
          fac=Lsup(dd_alpha)
          if norm_alpha > 0: fac*=1./norm_alpha
          error=fac*0.5*dt**2
          print "\testimated local error for time step size %e is %e"%(dt,error)
          dt_new=sqrt(2./3.*ODE_TOL*2/fac)
          dt_new=min(max(dt_new,dt/5),dt*5,DT_MAX) # avid rapit changes
          print "\tINFO: new time step size %e"%dt_new
    dt, dt_old=dt_new, dt
