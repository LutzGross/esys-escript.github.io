
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

# You can shorten the execution time by reducing variable tend from 60 to 0.5

# Importing all the necessary modules required.
from esys.escript import *
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
from numarray import identity,zeros,ones
import os
from phones import *

########################################################
# subroutine: wavesolver2d
# Can solve a generic 2D wave propagation problem with a
# point source in a homogeneous medium.
# Arguments:
#	domain  : domain to solve over
#	h       : time step
#	tend    : end time
#	lam, mu : lames constants for domain
#	rho	: density of domain
#	U0	: magnitude of source
#	xc	: source location in domain (Vector)
#	savepath: where to output the data files
########################################################
def wavesolver2d(domain,h,tend,lam,mu,rho,U0,xc,savepath):
   x=domain.getX()
   # ... open new PDE ...
   mypde=LinearPDE(domain)
   #mypde.setSolverMethod(LinearPDE.LUMPING)
   mypde.setSymmetryOn()
   kmat = kronecker(domain)
   mypde.setValue(D=kmat*rho)

   # define small radius around point xc
   # Lsup(x) returns the maximum value of the argument x
   src_radius = 50#2*Lsup(domain.getSize())
   print "src_radius = ",src_radius

   dunit=numarray.array([1.,0.]) # defines direction of point source

   
   # ... set initial values ....
   n=0
   # initial value of displacement at point source is constant (U0=0.01)
   # for first two time steps
   u=U0*(cos(length(x-xc)*3.1415/src_radius)+1)*whereNegative(length(x-xc)-src_radius)*dunit
   u_m1=u
   t=0

   u_pot = cbphones(domain,u,[[125,250],[250,250],[250,375]],2)
   print u_pot[0,1]
   u_pc_x = u_pot[0,0]
   u_pc_y = u_pot[0,1]

   # open file to save displacement at point source
   u_pc_data=open(os.path.join(savepath,'U_pc.out'),'w')
   u_pc_data.write("%f %f %f\n"%(t,u_pc_x,u_pc_y))
 
     #NEW WAY
     #g=grad(u)
     #stress=lam*trace(g)*kmat+mu*(g+transpose(g))

   while t<tend:
     # ... get current stress ....
     
##OLD WAY
     g=grad(u)
     stress=lam*trace(g)*kmat+mu*(g+transpose(g))
     ### ... get new acceleration ....
     #mypde.setValue(X=-stress)          
     #a=mypde.getSolution()
     ### ... get new displacement ...
     #u_p1=2*u-u_m1+h*h*a
###NEW WAY
     mypde.setValue(X=-stress*(h*h))
     mypde.setValue(Y=(rho*2*u-rho*u_m1))
     u_p1 = mypde.getSolution()
     # ... shift displacements ....
     u_m1=u
     u=u_p1
     #stress = 
     t+=h
     n+=1
     print n,"-th time step t ",t
     u_pot = cbphones(domain,u,[[125.,250.],[250.,250.],[250.,375.]],2)

#     print "u at point charge=",u_pc
     
     u_pc_x = u_pot[0,0]
     u_pc_y = u_pot[0,1]           
     # save displacements at point source to file for t > 0
     u_pc_data.write("%f %f %f\n"%(t,u_pc_x,u_pc_y))
 
     # ... save current acceleration in units of gravity and displacements 
     #saveVTK(os.path.join(savepath,"usoln.%i.vtu"%n),acceleration=length(a)/9.81,
     #displacement = length(u), tensor = stress, Ux = u[0] )
     saveVTK(os.path.join(savepath,"tonysol.%i.vtu"%n),output1 = length(u),tensor=stress)

   u_pc_data.close()