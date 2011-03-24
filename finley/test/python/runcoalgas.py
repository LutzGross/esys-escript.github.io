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
Coal Seam gasL ECLIPSE test case
"""
__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript import unitsSI as U
from coalgas import *
import time
from esys.finley import Rectangle

L_X=168*U.km
L_Y=168*U.km
L_Z=10*U.m

N_X=21 
N_Y=21 


PERM_F_X = 100 * U.mDarcy
PERM_F_Y = 100 * U.mDarcy
PERM_F_Z = 1e-4 * U.mDarcy
PHI_F_0 = 0.01
P_F_0 = 69 * U.bar

# these object correspond to the ECLIPSE input files 
PVTW={ "p_ref" :   1000 * U.bar ,  
       "B_ref" :  0.997  ,
       "C" :  3.084E-06  /U.bar,
       "mu_ref" : 0.68673 * U.cPoise,
       "C_v" : 0/U.bar

     }
GRAVITY = { "water" : 1.0, 
            "gas" : .553 }
ROCK = { "p_ref" :   1000 * U.bar ,
         "C" : 3.3E-004 * 1./U.bar }
     
LANGMUIR = [
[ 0	* U.bar , 0.00000000 ],
[ 100	* U.bar , 0.00213886 ],
[ 200	* U.bar , 0.00383259 ],
[ 300	* U.bar , 0.00520706 ],
[ 400	* U.bar , 0.00634474 ],
[ 500	* U.bar , 0.00730199 ],
[ 600	* U.bar , 0.00811857 ],
[ 700	* U.bar , 0.00882336 ],
[ 800	* U.bar , 0.00943786 ],
[ 900	* U.bar , 0.00997836 ],
[ 1000	* U.bar , 0.01045748 ],
[ 1200	* U.bar , 0.01126912 ] ]

PVDG = [
[ 14.70 * U.bar ,200.3800 , 0.012025 * U.cPoise ] ,
[ 20.00 * U.bar ,146.0600 , 0.012030 * U.cPoise ] ,
[ 25.00 * U.bar ,116.1461 , 0.012034 * U.cPoise ] ,
[ 30.00 * U.bar ,96.3132 , 0.012038 * U.cPoise ] ,
[ 35.00 * U.bar ,82.2113 , 0.012043 * U.cPoise ] ,
[ 49.33 * U.bar ,57.7891 , 0.012055 * U.cPoise ] ,
[ 59.00 * U.bar ,48.0866 , 0.012064 * U.cPoise ] ,
[ 69.00 * U.bar ,40.9441 , 0.012073 * U.cPoise ] ,
[ 75.00 * U.bar ,37.5839 , 0.012078 * U.cPoise ] ,
[ 83.00 * U.bar ,33.8685 , 0.012085 * U.cPoise ] ,
[ 90.00 * U.bar ,31.1661 , 0.012092 * U.cPoise ] ,
[ 95.00 * U.bar ,29.4827 , 0.012097 * U.cPoise ] ,
[ 100.00 * U.bar ,27.9698 , 0.012101 * U.cPoise ] ,
[ 105.00 * U.bar ,26.6028 , 0.012106 * U.cPoise ] ,
[ 118.60 * U.bar ,23.4749 , 0.012119 * U.cPoise ] ,
[ 120.00 * U.bar ,23.1937 , 0.012120 * U.cPoise ] ,
[ 140.00 * U.bar ,19.7977 , 0.012140 * U.cPoise ] ,
[ 153.23 * U.bar ,18.0443 , 0.012153 * U.cPoise ] ,
[ 160.00 * U.bar ,17.2607 , 0.012159 * U.cPoise ] ,
[ 170.00 * U.bar ,16.2188 , 0.012169 * U.cPoise ] ,
[ 187.86 * U.bar ,14.6373 , 0.012188 * U.cPoise ] ,
[ 222.49 * U.bar ,12.3027 , 0.012224 * U.cPoise ] ,
[ 257.13 * U.bar ,10.6038 , 0.012262 * U.cPoise ] ,
[ 291.76 * U.bar ,9.3134 , 0.012301 * U.cPoise ] ,
[ 326.39 * U.bar ,8.3001 , 0.012341 * U.cPoise ] ,
[ 361.02 * U.bar ,7.4835 , 0.012383 * U.cPoise ] ,
[ 395.66 * U.bar ,6.8114 , 0.012425 * U.cPoise ] ,
[ 430.29 * U.bar ,6.2491 , 0.012470 * U.cPoise ] ,
[ 464.92 * U.bar ,5.7715 , 0.012515 * U.cPoise ] ,
[ 499.55 * U.bar ,5.3610 , 0.012562 * U.cPoise ] ,
[ 534.19 * U.bar ,5.0043 , 0.012610 * U.cPoise ] ,
[ 568.82 * U.bar ,4.6917 , 0.012659 * U.cPoise ] ,
[ 603.45 * U.bar ,4.4154 , 0.012710 * U.cPoise ] ,
[ 638.08 * U.bar ,4.1695 , 0.012762 * U.cPoise ] ,
[ 672.72 * U.bar ,3.9491 , 0.012815 * U.cPoise ] ,
[ 707.35 * U.bar ,3.7507 , 0.012869 * U.cPoise ] ,
[ 741.98 * U.bar ,3.5711 , 0.012925 * U.cPoise ] ,
[ 776.61 * U.bar ,3.4076 , 0.012982 * U.cPoise ] ,
[ 811.25 * U.bar ,3.2583 , 0.013041 * U.cPoise ] ,
[ 845.88 * U.bar ,3.1214 , 0.013100 * U.cPoise ] ,
[ 880.51 * U.bar ,2.9953 , 0.013161 * U.cPoise ] ,
[ 915.14 * U.bar ,2.8790 , 0.013223 * U.cPoise ] ,
[ 949.78 * U.bar ,2.7712 , 0.013287 * U.cPoise ] ,
[ 984.41 * U.bar ,2.6711 , 0.013352 * U.cPoise ] ,
[ 1019.00 * U.bar ,2.5781 , 0.013418 * U.cPoise ] ,
[ 1053.70 * U.bar ,2.4909 , 0.013486 * U.cPoise ] ,
[ 1088.30 * U.bar ,2.4096 , 0.013554 * U.cPoise ] ,
[ 1122.90 * U.bar ,2.3334 , 0.013624 * U.cPoise ] ,
[ 1157.60 * U.bar ,2.2616 , 0.013696 * U.cPoise ] ,
[ 1192.20 * U.bar ,2.1942 , 0.013768 * U.cPoise ] ,
[ 1226.80 * U.bar ,2.1307 , 0.013842 * U.cPoise ] ,
[ 1261.50 * U.bar ,2.0705 , 0.013917 * U.cPoise ] ,
[ 1296.10 * U.bar ,2.0138 , 0.013994 * U.cPoise ] ,
[ 1330.70 * U.bar ,1.9600 , 0.014072 * U.cPoise ] ,
[ 1365.40 * U.bar ,1.9089 , 0.014151  * U.cPoise ] ]  


SGFN   = [  
[ 0  , 0  , 0 ],
[ 0.05   , 0  , 0 ],
[ 0.1333  , 0.00610   , 0 ],
[ 0.2167  , 0.02990   , 0 ],
[ 0.3  , 0.0759   , 0 ],
[ 0.3833  , 0.1471     , 0 ],
[ 0.46667  , 0.2458     , 0 ],
[ 0.55  , 0.3739     , 0 ],
[ 0.6333  , 0.53300    , 0 ],
[ 0.7167  , 0.7246     , 0 ],
[ 0.8  , 0.95       , 0 ] ]
SWFN = [
[ 0.20000  , 0.00000, 0 ],
[ 0.28330  , 0.03280, 0 ],
[ 0.36670  , 0.09270, 0 ],
[ 0.45000  , 0.17030, 0 ],
[ 0.53330  , 0.26220, 0 ],
[ 0.61670  , 0.36650, 0 ],
[ 0.70000  , 0.48170, 0 ], 
[ 0.78330  , 0.60710, 0 ],
[ 0.86670  , 0.74170, 0 ],
[ 0.95000  , 0.88500, 0 ],
[ 1.00000  , 1.00000, 0 ] ]


wellspecs = {
  'P1' : { "X0" : [106, 106], 
           "r"  : 0.253*U.m,
           "s"  : 0,
           "Q"  : 2000*U.Barrel/U.day*GRAVITY["water"],
           "BHP" : 75*U.psi,
           "schedule" : [0.*U.yr, 4*U.yr]
         }
}

CELL_X=L_X/N_X
CELL_Y=L_Y/N_Y
# print input
print("<%s> Execution started."%time.asctime())

print "length x-direction = %f km"%(L_X/U.km)
print "number of cells in x direction = %d"%N_X
print "cell size in x direction = %f m"%(CELL_X/U.m)
print "length y-direction = %f km"%(L_Y/U.km)
print "number of cells in y direction = %d"%N_Y
print "cell size in y direction = %f m"%(CELL_Y/U.m)
print "fracture permeability in x direction= %f mD"%(PERM_F_X/(U.mDarcy))
print "fracture permeability in y direction= %f mD"%(PERM_F_Y/(U.mDarcy))
print "fracture permeability in z direction= %f mD"%(PERM_F_Z/(U.mDarcy))
print "initial porosity in fractured rock= %f"%PHI_F_0

domain=Rectangle(N_X, N_Y, l0=L_X, l1=L_Y)
print("<%s> Mesh set up completed."%time.asctime())


model = PorosityOneHalfModel(domain, 
                             phi_f=Porosity(phi_0=PHI_F_0, p_0=P_F_0, p_ref=ROCK["p_ref"], C = ROCK["C"]),
                             L_g=InterpolationTable([ l[0] for l in LANGMUIR ], [ l[1] for l in LANGMUIR ] ),
			     perm_f_0=PERM_F_X, 
			     perm_f_1=PERM_F_Y, 
			     perm_f_2=PERM_F_Z,
			     k_w =InterpolationTable([ l[0] for l in SWFN ], [ l[1] for l in SWFN ] ),  
			     k_g= InterpolationTable([ l[0] for l in SGFN ], [ l[1] for l in SGFN ] ),  
			     mu_w = WaterViscosity(mu_ref = PVTW["mu_ref"], p_ref=PVTW["p_ref"], C=PVTW["C_v"]),      
			     mu_g = InterpolationTable([ l[0] for l in PVDG ], [ l[2] for l in PVDG ] ),
			     rho_w = WaterDensity(B_ref=PVTW["B_ref"], p_ref = PVTW["p_ref"], C=PVTW["C"], gravity=GRAVITY["water"]), 
			     rho_g=GasDensity( p = [ l[0] for l in PVDG ], B = [ l[1] for l in PVDG ], gravity=GRAVITY["gas"]), 
			     wells=[ VerticalPeacemanWell(i,BHP=wellspecs[i]["BHP"], 
							    Q=wellspecs[i]["Q"], 
							    r=wellspecs[i]["r"], 
							    X0=[ (wellspecs[i]["X0"][0]+0.5)*CELL_X,  (wellspecs[i]["X0"][1]+0.5)*CELL_Y],
							    D=[CELL_X, CELL_Y, L_Z], 
							    perm=[PERM_F_X, PERM_F_Y, PERM_F_Z], 
							    schedule=wellspecs[i]["schedule"], 
							    s=wellspecs[i]["s"]) for i in wellspecs] 
			     )
			     
model.setInitialState(p=P_F_0, S_fg=0,  C_mg=None)


1/0


#=======================================================

from esys.escript.linearPDEs import LinearPDE
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
DT=1.*U.sec

CHI=1.0   # =0 = explicit Euler, 0.5 = Crank-Nicholson, 1. = implicit Euler

G=-9.81*U.m/U.sec**2
S_W_MIN = 0.01
S_W_MAX = 0.75

M=3./2.                   *0+1
GAMMA=1.e-4*(1./U.Pa)

V_L=0.015*U.m**3/U.kg   *0
P_L=6.109*U.Mega*U.Pa
P_A=101.325*U.Kilo*U.Pa
RHO_AG=0.717*U.kg/U.m**3
P_0G=P_A
RHO_0G=RHO_AG
RHO_W = 1000.*U.kg/U.m**3 
RHO_COAL=1250*U.kg/U.m**3
ETA_W=1.e-3*U.kg/U.m/U.sec
ETA_G=1.1e-5*U.kg/U.m/U.sec
LAM=1./3.
PHI_0=0.00804
EPS_L=0.02295
K_S=8425*U.Mega*U.Pa
ALPHA=0.005*1./U.Pa            *0
EPS_V=0
K_0=3.7996e-17*U.m**3
D_G=20e-10*0


TOL=1.e-4
OUT_DIR="ReservoirResults"


class Porosity(object):
    def __init__(self, phi0=1., alpha=0., Ks=1., eps_v=0, p0=0., eps_L=0, p_L=0):
          self.phi0=phi0
          self.alpha=alpha
          self.eps_v=eps_v
          self.p_0=p0
          self.p_L=p_L
          self.eps_L = eps_L
          self.Ks = Ks

    def getValue(self,p=0):
         return self.phi0+ self.alpha*(self.eps_v + (p-self.p_0)/self.Ks + self.eps_L * self.p_L * (self.p_0-p)/((self.p_L + self.p_0) * (p + self.p_L)))

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
    def __init__(self):
         pass

    def __call__(self, p=0):
         return self.getValue(p)

    def getValue(self,p=0):
           return 0.*p
    def getDp(self, p=0):
         return 0.*p

class GasStorageCapacity(object):
    def __init__(self):
        raise NotImplementedError

    def getValue(self,p=0):
         raise NotImplementedError
    def getDp(self, p=0):
         raise NotImplementedError
    def __call__(self, p=0):
         return self.getValue(p)



class LangmuirIsotherm(GasStorageCapacity):
    def __init__(self, V_L=0, p_L=0., rho_a=1. ):
             self.V_L=V_L
             self.p_L=p_L
             self.rho_a=rho_a

    def getValue(self,p=0):
             return (self.V_L*self.rho_a)*p/(self.p_L +p)

    def getDp(self, p=0):
             return (self.V_L*self.rho_a)*p/(self.p_L +p)**2

  
class WaterVolumeFraction(object):
    """
    class defining water volume fraction S_w as function of suction pressure p
    """
    def __init__(self):
         raise NotImplementedError
    def __call__(self, p=0):
         return self.getValue(p)
    def getValue(self,p=0):
         """
         returns water volume fraction for given suction pressure p
         """
         raise NotImplementedError
    def getDp(self, p=0):
         """
         returns derivative of water volume fraction with respect to suction pressure for given suction pressure p
         """
         raise NotImplementedError

class VanGenuchten(WaterVolumeFraction):
    def __init__(self, S_min=0., S_max=1., m =1, gamma=0):
               self.S_min = S_min
               self.S_max = S_max
               self.m = m
               self.gamma = gamma
    def getEffFraction(self,s):
              return (s-self.S_min)/(self.S_max-self.S_min)
    def getFraction(self,s_eff):
              return self.S_min+s_eff*(self.S_max-self.S_min)

    def getValue(self,p=0):
             s_hat = (1+(abs(p) * self.gamma)**self.m)**(1./self.m-1)
             return self.getFraction(s_hat)
    def getDp(self, p=0):
         return (self.S_max-self.S_min) * (self.m-1)/self.m*(1+(self.gamma * abs(p))**self.m)**(1./self.m) * self.gamma**self.m * abs(p)**(self.m-1)*sign(p)

if DIM==2:
    dom=Rectangle(NE_L,NE_H,l0=L,l1=H,order=1,optimize=True)
else:
    dom=Brick(NE_L,NE_L,NE_H,l0=L,l1=L,l2=H,order=1,optimize=True)



S_w=VanGenuchten(S_min=S_W_MIN, S_max=S_W_MAX, m =M, gamma=GAMMA)
C_c=LangmuirIsotherm(V_L=V_L, p_L=P_L, rho_a=RHO_AG )

Rho_g=Density(p_a = P_A, rho_a= RHO_AG, p0 = 2*P_A, rho0=RHO_0G)
Phi=Porosity(phi0=PHI_0, alpha=ALPHA, Ks=K_S, eps_v=EPS_V, p0=P_A, eps_L=EPS_L, p_L=P_L)

pde=LinearPDE(dom,numEquations=3, numSolutions=3)
pde.getSolverOptions().setVerbosityOn()
pde.getSolverOptions().setSolverMethod(pde.getSolverOptions().DIRECT)

t=0
n=0
dt=DT
mkDir(OUT_DIR)

z=dom.getX()[DIM-1]
z_top=sup(z)
z_bot=inf(z)


suction=(RHO_W-RHO_AG)*G*(z-z_top)
p_g=RHO_AG*G*(z-z_top)+P_A
m_g=( 1-S_w(suction)) * ( Rho_g(p_g) * Phi(p_g)  - RHO_COAL * C_c(p_g) )

q=pde.createCoefficient("q")
r=pde.createCoefficient("r")
# fix pressures and movable mass on the surface
q[0]=whereZero(z-z_top)+whereZero(z-z_bot)
q[1]=whereZero(z-z_top)+whereZero(z-z_bot)
q[2]=whereZero(z-z_top)+whereZero(z-z_bot)
r[0]=suction
r[1]=m_g
r[2]=p_g
pde.setValue(q=q, r=r)

norm_p_g=Lsup(p_g)
norm_suction=Lsup(suction)

print "   m_g %e"%(Lsup(m_g))
print "   p_g %e"%(norm_p_g)
print "   suction %e"%(norm_suction)

while  t< T_END:
    suction_old=suction
    p_g_old=p_g
    m_g_old=m_g
    norm_dp_g = norm_p_g
    norm_dsuction = norm_suction

    print "Time step %d at %e sec"%(n+1,t+dt)
    niter=0

    while norm_dp_g > norm_p_g * TOL or norm_dsuction > norm_suction * TOL:
       print "   Iteration step %d:"%(niter,)

       s_w=S_w(suction)
       DS_wDp=S_w.getDp(suction)

       rho_g=Rho_g(p_g)
       DRho_gDp=Rho_g.getDp(p_g)

       c = C_c(p_g)
       DC_cDp=C_c.getDp(p_g)

       phi=Phi(p_g)
       DPhiDp=Phi.getDp(p_g)

       print "   porosity range phi:",inf(phi), sup(phi)
       print "   water volume fraction :",inf(s_w), sup(s_w)

       D=pde.createCoefficient("D")
       A=pde.createCoefficient("A")
       X=pde.createCoefficient("X")
       Y=pde.createCoefficient("Y")


       # p_g equation:
       Hp    =phi*DRho_gDp + rho_g*DPhiDp - RHO_COAL*DC_cDp
       Homega=DS_wDp*(rho_g*phi - RHO_COAL*c)
       D[2,0]=-Homega
       D[2,1]=-1.
       D[2,2]= Hp
       Y[2]=Hp * p_g_old - m_g_old - Homega*suction_old

       # revise :::
       k=K_0*(phi/PHI_0)**3
       # s_w_eff=S_w.getEffFraction(s_w)
       s_w_eff=s_w
       k_hat=sqrt(s_w_eff) * (1-(1-s_w_eff**(1./LAM))**LAM)**2
       k_w=k*k_hat
       k_g=k*(1-k_hat)

       K00=CHI*dt*k_w/ETA_W
       K11=CHI*dt*rho_g*D_G
       K12=CHI*dt*rho_g*k_g/ETA_G
 
       # suction equation:
       Y[0]=phi*DS_wDp*p_g_old+s_w*DPhiDp*suction_old
       X[0,:] = - dt * k_w/ETA_W * ( (1-CHI) * grad(suction+p_g) - G * RHO_W * kronecker(DIM)[DIM-1] )
       print Y[0]
       print X[0,:]
       D[0,0]=phi*DS_wDp
       D[0,2]=s_w*DPhiDp
       for i in range(DIM):
          A[0,i,0,i]=K00
          A[0,i,2,i]=K00


       # p_g equation:
       Y[1]=m_g_old
       X[1,:] = - dt * rho_g * ((1-CHI) * ( D_G * grad(m_g_old) + rho_g * k_g/ETA_G * ( grad(p_g_old) - G * rho_g * kronecker(DIM)[DIM-1]) ))
       print X[1,:]
       D[1,1]=1.
       for i in range(DIM):
          A[1,i,1,i]=K11
          A[1,i,2,i]=K12


       pde.setValue(A=A, D=D, X=X, Y=Y)
       u=pde.getSolution()
       m_g, m_g2=u[1], m_g
       p_g, p_g2=u[2], p_g
       suction, suction2=u[0], suction

       norm_p_g=Lsup(p_g)
       norm_suction=Lsup(suction)
       norm_dp_g = Lsup(p_g-p_g2)
       norm_dsuction = Lsup(suction- suction2)
       print inf(u[0])
       print sup(u[0])

       print "   m_g correction %e/%e"%(Lsup(m_g-m_g2), Lsup(m_g))
       print "   p_g correction %e/%e"%(norm_dp_g, norm_p_g)
       print "   suction correction %e/%e"%(norm_dsuction, norm_suction)

       niter+=1
       if niter >1: 1/0

    
    print "Time step %d completed."%n
    n+=1
    t+=dt

1/0

      

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
