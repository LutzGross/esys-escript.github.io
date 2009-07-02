
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Brick
from numpy import identity,zeros,ones

ne=32          # number of cells in x_0 and x_1 directions
width=10000.  # length in x_0 and x_1 directions
lam=3.462e9
mu=3.462e9
rho=1154.
tend=20. # to ran a full simulation change tend to 60.
alpha=0.7
t0=3.

U0=1. # maximum displacement
mkDir("data") # create directory data if it does not exist already.

def wavePropagation(domain,h,tend,lam,mu,rho, xc, src_radius, U0):
   x=domain.getX()
   # ... open new PDE ...
   mypde=LinearPDE(domain)
   mypde.getSolverOptions().setSolverMethod(mypde.getSolverOptions().LUMPING)
   kronecker=identity(mypde.getDim())

   dunit=numpy.array([1.,0.,0.]) # defines direction of point source

   mypde.setValue(D=kronecker*rho, q=whereNegative(length(x-xc)-src_radius)*dunit)
   # ... set initial values ....
   n=0
   # for first two time steps
   u=Vector(0.,Solution(domain))
   u_last=Vector(0.,Solution(domain))
   t=0

   # define the location of the point source 
   L=Locator(domain,xc)
   # find potential at point source
   u_pc=L.getValue(u)
   print "u at point charge=",u_pc
   # open file to save displacement at point source
   u_pc_data=FileWriter('./data/U_pc.out')
   u_pc_data.write("%f %f %f %f\n"%(t,u_pc[0],u_pc[1],u_pc[2]))
 
   while t<tend:
     t+=h
     # ... get current stress ....
     g=grad(u)
     stress=lam*trace(g)*kronecker+mu*(g+transpose(g))
     # ... get new acceleration ....
     amplitude=U0*(4*(t-t0)**3/alpha**3-6*(t-t0)/alpha)*sqrt(2.)/alpha**2*exp(1./2.-(t-t0)**2/alpha**2)
     mypde.setValue(X=-stress, r=dunit*amplitude)
     a=mypde.getSolution()
     # ... get new displacement ...
     u_new=2*u-u_last+h**2*a
     # ... shift displacements ....
     u_last=u
     u=u_new
     n+=1
     print n,"-th time step t ",t
     u_pc=L.getValue(u)
     print "u at point charge=",u_pc
     # save displacements at point source to file for t > 0
     u_pc_data.write("%f %f %f %f\n"%(t,u_pc[0],u_pc[1],u_pc[2]))
 
     # ... save current acceleration in units of gravity and displacements 
     if n==1 or n%10==0: saveVTK("./data/usoln.%i.vtu"%(n/10),acceleration=length(a)/9.81,
     displacement = length(u), tensor = stress, Ux = u[0] )

   u_pc_data.close()
  
mydomain=Brick(ne,ne,10,l0=width,l1=width,l2=10.*width/ne)
h=inf(1./5.)*inf(sqrt(rho/(lam+2*mu))*mydomain.getSize())
print "time step size = ",h
#  spherical source at middle of bottom face
xc=[width/2.,width/2.,0.]
# define small radius around point xc
src_radius = 0.03*width
print "src_radius = ",src_radius
wavePropagation(mydomain,h,tend,lam,mu,rho, xc, src_radius, U0)

