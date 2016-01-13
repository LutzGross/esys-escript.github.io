from esys.escript import *
from esys.escript.models import Mountains
from esys.finley import Brick,Rectangle
from math import pi, ceil

NE=16
DIM=3
H=1.
L=2*H
TOL=1.e-4
OMEGA=10
EPS=0.01
t=0
T_END=(2*pi)/OMEGA
n=0
if DIM==2:
  mydomain=Rectangle(int(ceil(L*NE/H)),NE,l0=L,l1=H,order=1, useFullElementOrder=True,optimize=True)
else:
  mydomain=Brick(int(ceil(L*NE/H)),int(ceil(L*NE/H)),NE,l0=L,l1=L,l2=H,order=1, useFullElementOrder=True,optimize=True)
x=mydomain.getX()
v = Vector(0.0, Solution(mydomain))
if DIM==2:
    a0=1
    n0=1
    n1=0.5
    a1=-(a0*n0)/n1
    v[0]=a0*sin(pi*n0*x[0])* cos(pi*n1*x[1])
    v[1]=a1*cos(pi*n0*x[0])* sin(pi*n1*x[1])
else:
    a0=1
    a1=1
    n0=2
    n1=2
    n2=0.5
    a2=-(a0*n0+a1*n1)/n2
    v[0]=a0*sin(pi*n0*x[0])* cos(pi*n1*x[1])* cos(pi*n2*x[2])
    v[1]=a1*cos(pi*n0*x[0])* sin(pi*n1*x[1])* cos(pi*n2*x[2])
    v[2]=a2*cos(pi*n0*x[0])* cos(pi*n1*x[1])* sin(pi*n2*x[2])


H_t=Scalar(0.0, Solution(mydomain))
mts=Mountains(mydomain,v,eps=EPS,z=H)
dt=0.
while t<T_END:
    print "STEP ", t
    u=v*cos(OMEGA*t)
    u,Z=mts.update(u=u,H_t=H_t,verbose=False)
    
    saveVTK("state.%d.vtu"%n,sol=Z)
    print "Integral(Z)=",integrate(Z),Lsup(u[DIM-1])
    n+=1
    H_t=Z
    dt=mts.getDt()
    t+=dt

