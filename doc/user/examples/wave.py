# $Id$
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Brick
ne=5           # number of cells in x_0-direction
depth=10000.   # length in x_0-direction
width=100000.  # length in x_1 and x_2 direction
lam=3.462e9
mu=3.462e9
rho=1154.
tau=10.
umax=2.
tend=60
h=1./5.*sqrt(rho/(lam+2*mu))*(depth/ne)
print "time step size = ",h

def s_tt(t): return umax/tau**2*(6*t/tau-9*(t/tau)**4)*exp(-(t/tau)**3)

def wavePropagation(domain,h,tend,lam,mu,rho,s_tt):
   x=domain.getX()
   # ... open new PDE ...
   mypde=LinearPDE(domain)
   mypde.setSolverMethod(mypde.LUMPING)
   mypde.setValue(D=kronecker(mypde.getDim())*rho, \
                  q=whereZero(x[0])*kronecker(mypde.getDim())[1,:])
   # ... set initial values ....
   n=0
   u=Vector(0,ContinuousFunction(domain))
   u_last=Vector(0,ContinuousFunction(domain))
   t=0
   while t<tend:
     # ... get current stress ....
     g=grad(u)
     stress=lam*trace(g)*kronecker(mypde.getDim())+mu*(g+transpose(g))
     # ... get new acceleration ....
     mypde.setValue(X=-stress,r=s_tt(t+h)*kronecker(mypde.getDim())[1,:])
     a=mypde.getSolution()
     # ... get new displacement ...
     u_new=2*u-u_last+h**2*a
     # ... shift displacements ....
     u_last=u
     u=u_new
     t+=h
     n+=1
     print n,"-th time step t ",t
     print "a=",inf(a),sup(a)
     print "u=",inf(u),sup(u)
     # ... save current acceleration in units of gravity
     if n%10==0: saveVTK("u.%i.xml"%(n/10),a=length(a)/9.81)

mydomain=Brick(ne,int(width/depth)*ne,int(width/depth)*ne,l0=depth,l1=width,l2=width)
wavePropagation(mydomain,h,tend,lam,mu,rho,s_tt)

