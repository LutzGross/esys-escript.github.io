from esys.escript import * 
from esys.escript.linearPDEs import LinearPDE, LinearSinglePDE
from esys.escript.nonlinearPDE import *
from esys.finley import Rectangle,Brick 
from math import pi
import os
import esys.escript.unitsSI as U 
#import numpy as np
#import scipy.optimize as so

#DIM=3
#NE=10

DIM=2
NE=80


class RockFeature(object):
     def __init__(self, lx, ly, lz, x=0, y=0, depth=0, rho=0):
         self.x=x
         self.y=y
         self.lx=lx
         self.ly=ly
         self.lz=lz
         self.depth=depth
         self.rho=rho
     def getMask(self, x):
       DIM=x.getDomain().getDim()
       m=whereNonPositive(x[DIM-1]+self.depth) * whereNonNegative(x[DIM-1]+(self.depth+self.lz)) \
         *whereNonNegative(x[0]-(self.x-self.lx/2)) * whereNonPositive(x[0]-(self.x+self.lx/2)) 
       if DIM>2:
        m*=whereNonNegative(x[1]-(self.y-self.ly/2)) * whereNonPositive(x[1]-(self.y+self.ly/2)) 
       return m	 
         
         
H=20*U.km
L=60*U.km
L0=0
L1=0
H_earth=10.*U.km


rho_rock=2300*U.kg/U.m**3
rho_air=0.
feastures = [ RockFeature(lx=10.*U.km, ly=10.*U.km, lz=1.*U.km, x= L/2-10.*U.km, y=L/2-10.*U.km, depth=1.*U.km, rho=rho_rock*0.1),
              RockFeature(lx=10.*U.km, ly=10.*U.km, lz=1.*U.km, x= L/2+10.*U.km, y=L/2+10.*U.km, depth=5.*U.km, rho=rho_rock*0.1) ]



              
# generate Domain:
NE_H=NE
NE_L=int((L/H)*NE+0.5)
if DIM==2:
   domain = Rectangle(NE_L,NE_H,l0=L,l1=H)
   x_cord=domain.getX()-[L0, H_earth]
else:   
   domain = Brick(NE_L,NE_L,NE_H,l0=L,l1=L,l2=H)
   x_cord=domain.getX()-[L0, L1, H_earth]
LL=max(L,H)
   
m_psi_ref=whereZero(x_cord[DIM-1]-inf(x_cord[DIM-1])) +  whereZero(x_cord[DIM-1]-sup(x_cord[DIM-1]))
for i in range(DIM-1):
    m_psi_ref= m_psi_ref + whereZero(x_cord[i]-inf(x_cord[i])) +  whereZero(x_cord[i]-sup(x_cord[i]))

# create test density:
rho_ref= 0
for f in feastures:
     m=f.getMask(x_cord)
     rho_ref = rho_ref * (1-m) + f.rho * m

#rho_ref=sup(x_cord[DIM-1])-x_cord[DIM-1] for testing!

# get the reference potential:
pde=LinearSinglePDE(domain)
pde.setValue(A=kronecker(domain), Y=4*pi*rho_ref, q=m_psi_ref)
pde.getSolverOptions().setVerbosityOn()
pde.setSymmetryOn()
#pde.getSolverOptions().setSolverMethod(pde.getSolverOptions().DIRECT) 
psi_ref=pde.getSolution()
del pde
d_obs=kronecker(DIM)[DIM-1]
g_hat=grad(psi_ref)[DIM-1]
beta=1/1000000.
#
#  where do we know the gravity:
#
x=Function(domain).getX()
dz=H/NE_H
H_earth=int(H_earth/dz)*dz # lock to grid
chi=whereNegative(abs(x[DIM-1]-(H_earth+dz/2))-dz/2)
# 
#   normalize  g_hat (data):
#
g0=Lsup(chi * g_hat)
if not g0 > 0: g0=1.
g_hat*=1./g0
print "Data normalization factor  = %e"%g0

#===============================================================================================================
#
#   use a linear PDE to solve the problem:
#
print "====== Use linear PDE ============================================="
pde=LinearPDE(domain, numEquations=3, numSolutions=3)
A=pde.createCoefficient("A")
X=pde.createCoefficient("X")
D=pde.createCoefficient("D")
q=pde.createCoefficient("q")


A[0,:,0,:]=kronecker(DIM) *beta 
A[1,:,1,:]=kronecker(DIM)
A[2,:,2,:]=kronecker(DIM)
A[2,:,1,:]=chi * outer(d_obs, d_obs)
D[0,2]=-4*pi/LL**2
D[1,0]=-4*pi/LL**2
X[2,:]= chi * g_hat * d_obs
#q[0]=wherePositive(domain.getX()[DIM-1]-H_earth)
#q[0]=whereZero(x_cord[DIM-1]-inf(x_cord[DIM-1])) +  whereZero(x_cord[DIM-1]-sup(x_cord[DIM-1]))
q[0]=wherePositive(domain.getX()[DIM-1]-H_earth) + whereZero(x_cord[DIM-1]-inf(x_cord[DIM-1])) +  whereZero(x_cord[DIM-1]-sup(x_cord[DIM-1]))
#q[0]=whereZero(x_cord[DIM-1]-sup(x_cord[DIM-1])) 
q[1]=m_psi_ref
q[2]=m_psi_ref

pde.setValue(A=A, D=D, X=X, q=q)
pde.getSolverOptions().setVerbosityOn()
pde.getSolverOptions().setTolerance(1e-8)
pde.getSolverOptions().setSolverMethod(pde.getSolverOptions().DIRECT) 

u=pde.getSolution()
rho_l, psi_l, =u[0], u[1]
print "rho =",rho_l
print "psi =",psi_l
print "lambda =",u[2]


#=========== This is the same with the variational class ===================================
print "====== Use variational problem  ============================================="
psi_s=Symbol("psi", (), dim=DIM)
rho_s=Symbol("rho", (), dim=DIM)
#g=Symbol("g", (), dim=DIM)

v=VariationalProblem(domain, u=psi_s,p=rho_s, debug=VariationalProblem.DEBUG3)
v.setValue( H = 0.5*chi*(grad(psi_s)[DIM-1]-g_hat)**2 + 0.5* beta * length(grad(rho_s))**2,
            X=grad(psi_s), Y=-4*pi*rho_s/LL**2,
            qp=q[0], q=m_psi_ref)
v.getNonlinearPDE().getLinearSolverOptions().setSolverMethod(pde.getSolverOptions().DIRECT) 
            
rho_v, psi_v, lag=v.getSolution(psi=0, rho=0)  # gamma=1 is the interesting case!
print "rho =",rho_v
print "rho_ref =",rho_ref*g0
print "psi =",psi_v
print "lambda =",lag

print "Differences PDE <-> Variational:"
print "rho =",Lsup(rho_v-rho_l)/Lsup(rho_l)
print "psi =",Lsup(psi_v-psi_l)/Lsup(psi_l)
print "Differences to Data :"
print "rho =",Lsup(rho_ref-rho_l*g0/LL**2)/Lsup(rho_ref)
print "psi =",Lsup(psi_ref-psi_l*g0)/Lsup(psi_ref)
print "g =",Lsup(chi*grad(psi_ref-psi_l*g0)[DIM-1])/Lsup(chi*grad(psi_ref)[DIM-1])


#saveVTK("u.vtu", rho_ref=rho_ref, psi_ref=psi_ref, g=grad(psi_ref)[DIM-1], rho=rho_v, psi=psi_v, chi=chi)
