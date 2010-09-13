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
DIM=2                           # spatial dimension
H=3.*U.m                     # height
L=2*H                           # length
NE_H=30                           # number of elements in H-direction. 
NE_L=int(ceil(L*NE_H/H))



C_V = 3e-11/U.Pa 
C_D = 5/U.sec
C_1 = 1e-12/U.sec
C_2 = 0.03
XI_0 = -0.8
LAME_0 = 46e9*U.Pa
MU_0 = 30e9*U.Pa
ALPHA_0 = 0
RHO = 2800*U.kg/U.m**3
G = 10*U.m/U.sec**2     *0
DT=1.*U.sec
VMAX=-0.01*U.m/U.sec

T_END=100000*U.yr                       # end time

VERBOSE=True
DT_VIS=T_END/100                # time distane between two visulaization files
DN_VIS=5                        # maximum counter increment between two visulaization files
VIS_DIR="results"               # name of the director for vis files

ODE_TOL=0.01
ODE_ITER_TOL=1.e-8
ODE_ITER_MAX=15

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
alpha=Scalar(0.,Function(dom))
u=Vector(0.,ContinuousFunction(dom))

pde=LinearPDESystem(dom)
pde.setSymmetryOn()
pde.getSolverOptions().setSolverMethod(pde.getSolverOptions().DIRECT)



fixed_v_mask=Vector(0,Solution(dom))
u0=Vector(0.,ContinuousFunction(dom))
for d in range(DIM):
    fixed_v_mask+=whereZero(x[d]-BBOX[d][0])*unitVector(d,DIM)
    if d == 0:
       fixed_v_mask+=whereZero(x[d]-BBOX[d][1])*unitVector(d,DIM)
       u0[d]=(x[d]-BBOX[d][0])/(BBOX[d][1]-BBOX[d][0])*VMAX

pde.setValue(Y=-G*RHO*kronecker(DIM)[DIM-1], q=fixed_v_mask, r=u0)
#
#  let the show begin:
#
k3=kronecker(DIM)
k3Xk3=outer(k3,k3)
alpha_old=alpha
dt_old=None
while t<T_END:

    print "start time step %d"%(n+1,)
    if n>1:
       print" eps :"
       print sqI2**2
       print 2*mu*(sqI2-gamma/(2*mu))*eps_e[1,1]-gamma*(sqI2-lame/gamma*I1)*sqI2

    I1=trace(eps_e)
    sqI2=length(eps_e)
    xi=safeDiv(I1,sqI2)
    i_xi=safeDiv(sqI2,I1)
    print "\txi = [ %e, %e]"%(inf(xi),sup(xi))
    
    # update damage parameter:
    m=wherePositive(xi-XI_0)
    a=sqI2**2*(xi-XI_0)*(m*C_D + (1-m)* C_1)
    b=(1-m)*(1./C_2)

    alpha, alpha_old, alpha_oold =solveODE(alpha, a,b, dt), alpha, alpha_old


    # step size for the next time step:

    print "\talpha = [ %e, %e]"%(inf(alpha),sup(alpha))
    if sup(alpha)  >1:
        raise ValueError,"damage parameter > 1 detected."

    gamma=alpha*GAMMA_M
    lame=LAME_0
    mu=MU_0*(1-alpha)

    lame_eff=lame-gamma*i_xi
    mu_eff=mu-gamma*xi
    print "\tmu_eff = [ %e, %e]"%(inf(mu_eff),sup(mu_eff))
    print "\tlame_eff = [ %e, %e]"%(inf(lame_eff),sup(lame_eff))
    if inf(mu_eff) <= 0: 
        raise ValueError,"Material failed due to non-positive mu_eff."

    i_eta = clip(2*C_V*(alpha-alpha_old)/dt, minval=0.)

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
          dt_new=min(max(dt_new,dt/5),dt*5) # avid rapit changes
          print "\tINFO: new time step size %e"%dt_new
    dt, dt_old=dt_new, dt
