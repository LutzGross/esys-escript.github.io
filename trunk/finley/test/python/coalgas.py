#######################################################
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
Gas in Coal Seam (fully coupled version)
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
from esys.weipa import saveVTK
from math import pi, ceil
import sys
import time

# ======================= Default Values ==================================================
H=2.*U.m                     # height
L=H                          # length
NE_H=5                       # number of elements in H-direction. 
NE_L=int(ceil(L*NE_H/H))
DIM=2
T_END=40000.0*U.sec                       	# end time  at max strain = 0.5

CHI=0.5   # =0 = explicit Euler, 0.5 = Crank-Nicholson, 1. = implicit Euler
G=-9.81*U.m/U.sec**2


class Porosity(object):
    def __init__(self, phi0=1., alpha=0., Ks=1., eps_v=0, p0=0., eps_L=0, p_L=0):
          self.phi0=phi0
          self.alpha=alpha
          self.eps_v=eps_v
          self.p0=p0
          self.p_L=p_L
          self.eps_L = eps_L
          self.Ks = Ks

    def getValue(self,p=0):
         return self.phi0*+ self.alpha*(self.eps_v + (p-self.p0)/Ks + self.eps_L * self.p_L * (self.p0-p)/((self.p_L + self.p_0) * (p + self.p_L)))

    def __call__(self, p=0):
         return self.getValue(p)

    def getDp(self, p=0):
          return self.alpha * (1./self.Ks - self.eps_L * self.p_L/(p-self.p_L)**2 )

class Density(object):
    def __init__(self, p_a =0, rho_a=0, p0 =0, rho0=0):
             self.p_a=p_a
             self.rho_a=rho_a
             self.p0=p0
             self.rho0=rho0
    def getValue(self,p=0):
             return (self.rho_a - self.rho0)/(self.p_a - self.p0)*(p-self.p0) + self.rho0
    def __call__(self, p=0):
         return self.getValue(p)

    def getDp(self, p=0):
          return (self.rho_a - self.rho0)/(self.p_a - self.p0)


class GasStorageCapacity(object):
    def __init__(self, V_L=0, p_L=0.):
             self.V_L=V_L
             self.p_L=p_L

    def getValue(self,p=0):
             return self.V_L*p/(self.p_L +p)
    def __call__(self, p=0):
         return self.getValue(p)

    def getDp(self, p=0):
         return self.V_L*p/(self.p_L +p)**2
  
class WaterVolumeFraction(object):
    def __init__(self, S_res=0.1, m =3, p_d=1.)
               self.S_res = S_res
               self.m = m
               self.p_d = p_d

    def getValue(self,p=0):
             s_hat = (1+(abs(p)/self.p_d)**self.m)**(1./self.m-1)
             return self.S_res+s_hat*(1-self.S_res)

    def __call__(self, p=0):
         return self.getValue(p)

    def getDp(self, p=0):
         tmp = (abs(p)/self.p_d)**self.m
         return (1-self.S_res) * (self.m-1)/self.m*((1+tmp)**(1./self.m) * tmp*1/p

if DIM==2:
    dom=Rectangle(NE_L,NE_H,l0=L,l1=H,order=1,optimize=True)
else:
    dom=Brick(NE_L,NE_L,NE_H,l0=L,l1=L,l2=H,order=1,optimize=True)



Rho_g=Density(p_a = PA, rho_a= RHO_AG, p0 =0, rho0=0)
S_w=WaterVolumeFraction(S_res=S_RES, m =M, p_d=P_D)
C_c=GasStorageCapacity(V_L=V_L, p_L=P_L)


pde=LinearPDE(dom,numEquations=3)

t=0
n=0

z0=sup(dom.getX()[DIM-1])
p_w=RHO_AW*G*(z-z0)+P_A
p_g=RHO_AG*G*(z-z0)+P_A
omega=p_w-p_g
m_g=(1-s_w) * ( rho_g * phi  + RHO_AG * RHO_COAL * c)


while :
    s_w=S_w(omega)
    rho_g=Rho_g(p_g)
    c = C_c(p_g)


    D=pde.getNewCoefficient("D")
    A=pde.getNewCoefficient("A")
    X=pde.getNewCoefficient("X")
    Y=pde.getNewCoefficient("Y")

    D[0,0]=phi*DS_w*Dp
    D[0,2]=s_w*DPhiDp

    D[1,1]=1.

    Hp    =phi*DRho_gDp+rho_g*DPhiDp+RHO_AG*RHO_COAL*DC_cDp
    Homega=Ds_wDo*(rho_g*phi+RHO_AG*RHO_COAL*c)
    D[2,0]=-Homega
    D[2,1]=-1.
    D[2,2]= Hp

 


    K00=chi*dt*k_w/ETA_W
    K11=chi*dt*rho_g*D_G
    K12=chi*dt*rho_g*k_g/ETA_G
    for i in range(DIM):
       A[0,i,0,i]=K00
       A[0,i,2,i]=K00
       A[1,i,1,i]=K11
       A[1,i,2,i]=K12


   Y[0]=phi*DS_w*Dp*p_g_old+s_w*DPhiDp*omega_old
   Y[1]=m_g_old
   Y[2]=h_p * p_g_old - m_g_old - Homega*omega


   X[0,:] = dt * k_w/ETA_W * ( (1-chi) * (grad(omega_old) + grad(omega_g_old) ) - g * RHO_W * kronecker(DIM)[DIM-1] )
   X[1,:] = dt * rho_g ((1-chi) * ( D_G * grad(m_g_old) + rho_g * k_g/eta_g * grad(p_g_old) - g * rho_g * k_g/eta_g * kronecker(DIM)[DIM-1]) )

#===============================
def meanValue(arg):
   return integrate(arg,Function(arg.getDomain()))/(H*L)

#Boundary conditions: 
#	axial loading: they applied a stress inversely proportional to the acoustic emission rate. We could have the axial forcing a stress or velocity inversely proportional to dalpha/dt (only when it is positive, and with the applied forcing rate going to zero when damage accumulation rate goes to a value we can determine in a test run with constant forcing). If this is to challenging or time consuming we could have a constant axial strain rate with very short time steps (at least when alpha increases above 0.3).

#Variables calculated and written to an output file:
#	time
#	differential stress (S_33-S_11)
#	deviatoric stress (S_33 - p)
#	Axial and transverse strain
#	damage and damage rate

	

#			    			# material parameters
C_V = 2.e-5/(U.Mega*U.Pa)  
C_D = 3/U.sec   			
C_1 = 1e-10/U.sec   			#1e-12/U.sec  
C_2 = 0.02				#0.03
XI_0 = -0.56			
LAME_0 = 29e9*U.Pa
MU_0 = 19e9*U.Pa
ALPHA_0 = 0.
ALPHA_pert = 0.6
RHO = 2800*U.kg/U.m**3
G = 10*U.m/U.sec**2     *0
#			    		# BC parameters
SIGMA_N=+50.*U.Mega*U.Pa		#confining -50 ,  +50 gives tension -> high damage levels
DIM=2                          
VMAX=-1.*U.m/U.sec/800000.		#*2
DT_MAX=50000.*U.sec
DT=400.*U.sec/10000.
#T_END=30*DT
if DIM==2:
  xc=[L/2,H/2]
else:
  xc=[L/2,L/2,H/2]
WWW=min(H,L)*0.02		 	# for one-element with a_0=pert+ALPHA_0 should be ~element length (=L/NE_L)

VERBOSE=True
DT_VIS=T_END/100                	# time distane between two visulaization files
DN_VIS=20                        	# maximum counter increment between two visulaization files
VIS_DIR="damage-localization"        	# name of the directory for vis files    		 <<<<<<<<<<<<<<<<<<<
DIAGNOSE_FN="diagnose_damage.csv"

ODE_TOL=0.1
ODE_ITER_TOL=1.e-6
ODE_ITER_MAX=15	
DEPS_MAX=0.01
DU_ITER_TOL=1e-4
DU_ITER_MAX=20

diagnose=FileWriter(DIAGNOSE_FN,append=False)
#===================================
S=0.5*XI_0*((2.*MU_0+3.*LAME_0)/(3.-XI_0**2) + LAME_0)
GAMMA_M=S + sqrt(S**2+2.*MU_0*(2.*MU_0+3.*LAME_0)/(3.-XI_0**2))

class DivergenceError(Exception):
    pass

def solveODE(u0, a, b, dt):
    """
    solves du/dt=a*exp(b*u)   u(t=0)=u0 and return approximation at t=dt
		
    we use backwards Euler by solving     u-u0 = dt * a*exp(b*u)		alpha=solveODE(alpha_old, a,b, dt)
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
         if n>ODE_ITER_MAX:
            print "*** ODE iteration has given up after %s steps with correction %e"%(n, norm_du)
            raise DivergenceError,"ODE iteration failed after %s steps."%n
    print "     ODE iteration completed after %d steps with correction %e."%(n,norm_du)
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
#   	........set up domain 
#

BBOX=boundingBox(dom)
DIM=dom.getDim()
x=dom.getX()
#
#   		initial values:
#
sigma=Tensor(0.,Function(dom))
eps_e=Tensor(0.,Function(dom))
#
#   		nucleation site - alpha perturbation:

alpha=ALPHA_0+ALPHA_pert*exp(-length(Function(dom).getX()-xc)**2/WWW**2)
#MU_0 = MU_0 + MU_0/100*exp(-length(Function(dom).getX()-xc)**2/WWW**2)

#saveVTK("alpha0_test-cd_0-mu0_pert.vtu",alpha=alpha)  # map initial damage.vtu		<<<<<<<<<<<<<<<<<<  


# pde.getSolverOptions().setPreconditioner(pde.getSolverOptions().AMG)
pde.getSolverOptions().setTolerance(DU_ITER_TOL**2)



fixed_v_mask=Vector(0,Solution(dom))
v0=Vector(0.,ContinuousFunction(dom))

for d in range(DIM):
    fixed_v_mask+=whereZero(x[d]-BBOX[d][0])*unitVector(d,DIM)
    if d == DIM-1:
        fixed_v_mask+=whereZero(x[d]-BBOX[d][1])*unitVector(d,DIM)
        v0[d]=(x[d]-BBOX[d][0])/(BBOX[d][1]-BBOX[d][0])*VMAX

pde.setValue(Y=-G*RHO*kronecker(DIM)[DIM-1], q=fixed_v_mask)
du=Vector(0.,Solution(dom))
u=Vector(0.,Solution(dom))
norm_du=0.
deps=Tensor(0,Function(dom))
i_eta=0
#
#  	........let the show begin:
#
k3=kronecker(DIM)
k3Xk3=outer(k3,k3)
alpha_old=alpha
dt_old=None
diagnose.write("t, -e22, e11, s00-s22, alpha_max, alpha_av\n")

while t<T_END:

    print "===== start time step %d ===== "%(n+1,)

    eps_e_old = eps_e
    sigma_old = sigma 
    alpha_old, alpha_oold = alpha, alpha_old
    du_old=du
    converged=False
    #  start the iteration for deps on a time step: deps from the last time step is used as an initial guess:
    while not converged:
        iter_du=0
        if dt_old!=None: 
             du=du_old*(dt/dt_old)
        else:
             du=du_old
        norm_ddu=Lsup(du)
        norm_u = Lsup(u)
        deps=symmetric(grad(du))
        sigma=sigma_old
        # start iteration :
        try: 
            while norm_ddu > DU_ITER_TOL * norm_u or iter_du == 0 :
               print "  start iteration step %d at time step %d:"%(iter_du,n+1)

               eps_e = eps_e_old + deps-(dt/2)*i_eta*deviatoric(sigma)
       
               I1=trace(eps_e)
               sqI2=length(eps_e)
               xi=safeDiv(I1,sqI2)
               # ......update damage parameter:
               m=wherePositive(xi-XI_0)
               a=sqI2**2*(xi-XI_0)*(m*C_D + (1-m)* C_1)
               b=(1-m)*(1./C_2)

               alpha=solveODE(alpha_old, a,b, dt)
               if sup(alpha) >  1.:
                  print "*** damage parameter %e > 1"%(sup(alpha), )
                  raise DivergenceError,"damage parameter %e > 1"%(sup(alpha), )

               alpha_dot=(alpha-alpha_old)/dt
               i_eta = clip(2*C_V*alpha_dot,minval=0.)
    
               gamma=alpha*GAMMA_M
               lame=LAME_0
               mu=MU_0*(1-alpha)

               print "     alpha = [ %e, %e]"%(inf(alpha),sup(alpha))
               print "     xi = [ %e, %e]"%(inf(xi),sup(xi))
               print "     gamma = [ %e, %e]"%(inf(gamma),sup(gamma))

               sigma = lame * I1 * k3 + 2* mu * eps_e - gamma * ( sqI2 * k3 + xi * eps_e )

               pde.setValue(A = mu * ( swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3) ) + lame*k3Xk3)
               pde.setValue(X=-sigma, y=SIGMA_N*dom.getNormal(), r=dt*v0-du)
        

               ddu=pde.getSolution()
               deps+=symmetric(grad(ddu))
               du=du+ddu
               norm_ddu=Lsup(ddu)
               norm_u=Lsup(u+du)
               print "  displacement change update = %e of %e"%(norm_ddu, norm_u)
               iter_du+=1
               if iter_du >  DU_ITER_MAX:
                  print "*** displacement iteration has given up after %s steps with rel. correction %e"%(iter_du, norm_ddu/norm_u)
                  raise DivergenceError,"displacement iteration failed after %s steps."%iter_du
            converged=True
            print "displacement iteration converged after %d steps (rel. increment = %e)."%(iter_du, norm_ddu/norm_u)

        except DivergenceError:
            converged=False
            dt*=0.5
            print "*** iteration is resumed with new time step size = %e"%(dt,)
              
    u+=du
    n+=1
    t+=dt
    diagnose.write(("%e,"*6+"\n")%(t, meanValue(-symmetric(grad(u))[DIM-1,DIM-1]), 
                                      meanValue(symmetric(grad(u))[0,0]),
                                      meanValue(sigma[0,0]-sigma[DIM-1,DIM-1]), 
                                      sup(alpha),
                                      meanValue(alpha)))	#meanValue(alpha_dot)))
    print "time step %s (t=%e) completed."%(n,t)
    #
    #  ............visualization:
    #
    # dom.setX(dom.getX()+du)
    if t>=t_vis or n>n_vis:
      saveVTK(os.path.join(VIS_DIR,"state.%d.vtu"%counter_vis),u=u, alpha=alpha,  I1=trace(eps_e), eps_12=(symmetric(grad(u))[0,1]), xi=safeDiv(trace(eps_e),length(eps_e)), alpha_dot=alpha_dot, sqI2=sqI2)	# tau=sigma[0,1], I2=length(eps_e)**2
      print "visualization file %d for time step %e generated."%(counter_vis,t)
      counter_vis+=1
      t_vis+=DT_VIS
      n_vis+=DN_VIS
    # 
    #   ............control time step size:
    #
    # ss=sup(length(deps))
    # if ss>0:
    #    dt_new=DEPS_MAX/ss*dt
    #    print "  time step size to control strain increment %s."%(dt_new,)
    # else:
    dt_new=dt
    if dt_old != None:
          dd_alpha=2.*dt_old*(alpha-alpha_old)+(alpha_oold-alpha_old)*dt/(dt*dt_old*(dt_old+dt))
          norm_alpha=Lsup(alpha)
          fac=Lsup(dd_alpha)
          if norm_alpha > 0: fac*=1./norm_alpha
          if fac>0:
             error=fac*0.5*dt**2
             print "  estimated local error for time step size %e is %e"%(dt,error)
             dt_new=sqrt(ODE_TOL*2/fac)
             print "  new time step size to maintain error level = %e."%(dt_new,)
             if n<10:
                fac=10.
             else:
                fac=1.3
             dt_new=min(max(dt_new,dt/fac),dt*fac) # aviod rapid changes				
             print "new time step size %e"%dt_new
    dt, dt_old=dt_new, dt
print "#### Time integartion completed #########"
