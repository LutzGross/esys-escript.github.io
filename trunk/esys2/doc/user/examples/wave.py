# $Id$
from esys.escript import *
import esys.finley
import numarray

def wavePropagation(domain,dt,tend,lame_lambda,lame_mu,rho,xc,r,tau,umax):
   # ... get characteristic function of impact region:
   x=domain.getX()
   location=x[0].whereZero()*(length(x-xc)-r).whereNegative()
   # ... open new PDE ...
   myPDE=LinearPDE(mydomain)
   myPDE.setLumpingOn()
   myPDE.setValue(D=numarray.identity(myPDE.getDim())*rho,q=location*numarray.identity(myPDE.getDim())[:,0])
   # ... set the initial values :
   t=dt
   n=0
   u=0
   v=0
   a=0
   while t<tend:
     # ... up-date displacement ....  
     u=u+dt*v+dt**2/2*a
     # ... get current stress ....
     g=grad(u)
     stress=lame_lambda*trace(g)+lame_mu*(g+transpose(g))
     # ... get new acceleration ....
     myPDE.setValue(X=stress,q=impact_location, \
                r=umax/tau**2*(6*t/tau-9*(t/tau)^4)*exp(-(t/tau)^3)*numarray.identity(myPDE.getDim())[:,0])
     a_new=myPDE.getSolution()
     # ... update velocity ... 
     v=v+h/2*(a+a_new)
     # ... next time step ...
     a=a_new
     t+=dt
     n+=1
     # ... save current displacement:
     if n%10: u.saveDX("u.%i.dx"%n)

ne=6 
lame_lambda=3.462e9
lame_mu=3.462e9
rho=1154.
tau=2.
umax=15.
xc=[0,1000,1000]
r=1.
tend=10.
dt=1./5.*sqrt(rho/(lame_lambda+2*lame_mu))(20000./ne)
print "step size = ",dt
mydomain=finely.Brick(ne,10*ne,10*ne,l0=20000,l1=200000,l2=200000)
wavePropagation(mydomain,dt,tend,lame_lambda,lame_mu,rho,xc,r,tau,umax)

