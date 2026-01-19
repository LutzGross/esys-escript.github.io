##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import matplotlib
matplotlib.use('agg')    #For interactive use, you can comment out this line
#It's just here to make testing easier
import matplotlib.pyplot as plt
from numpy import zeros,ones
import numpy
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.pdetools import Locator
try:
    from esys.finley import Brick
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False
from esys.weipa import saveVTK

if not HAVE_FINLEY:
    print("Finley module not available")
else:
    ne=32          # number of cells in x_0 and x_1 directions
    width=10000.  # length in x_0 and x_1 directions
    lam=3.462e9
    mu=3.462e9
    rho=1154.
    tend=10. # to ran a full simulation change tend to 60.
    alpha=0.7
    t0=3.


    U0=1. # maximum displacement
    mkDir("output") # create directory output if it does not exist already.

    def wavePropagation(domain,h,tend,lam,mu,rho, xc, src_radius, U0):
       # lists to collect displacement at point source
       ts, u_pc0,u_pc1,u_pc2=[], [], [], []
       x=domain.getX()
       # ... open new PDE ...
       mypde=LinearPDE(domain)
       mypde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
       kron=kronecker(mypde.getDim())

       dunit=numpy.array([1.,0.,0.]) # defines direction of point source

       mypde.setValue(D=kron*rho, q=whereNegative(length(x-xc)-src_radius)*dunit)
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
       print("u at point charge = %s"%u_pc)
       ts.append(t); u_pc0.append(u_pc[0]), u_pc1.append(u_pc[1]), u_pc2.append(u_pc[2])
     
       while t<tend:
         t+=h
         # ... get current stress ....
         g=grad(u)
         stress=lam*trace(g)*kron+mu*(g+transpose(g))
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
         print("time step %d, t = %s"%(n,t))
         u_pc=L.getValue(u)
         print("u at point charge = %s"%u_pc)
         ts.append(t); u_pc0.append(u_pc[0]), u_pc1.append(u_pc[1]), u_pc2.append(u_pc[2])
     
         # ... save current acceleration in units of gravity and displacements 
         if n==1 or n%10==0: saveVTK("output/usoln.%i.vtu"%(n/10),acceleration=length(a)/9.81,
         displacement = length(u), tensor = stress, Ux = u[0] )
       return ts, u_pc0,u_pc1,u_pc2

    #
    # create domain:
    #
    mydomain=Brick(ne,ne,10,l0=width,l1=width,l2=10.*width/ne)
    #
    #  sety time step size:
    #
    h=inf(1./5.)*inf(sqrt(rho/(lam+2*mu))*mydomain.getSize())
    print("time step size = %s"%h)
    #
    #  spherical source at middle of bottom face
    #
    xc=[width/2.,width/2.,0.]
    # define small radius around point xc
    src_radius = 0.03*width
    print("src_radius = %s"%src_radius)
    #
    # run it
    #
    ts, u_pc0,u_pc1,u_pc2 = wavePropagation(mydomain,h,tend,lam,mu,rho, xc, src_radius, U0)
    #
    # create a plot:
    #
    if getMPIRankWorld() == 0:
        plt.title("Displacement at Point Source")
        plt.plot(ts, u_pc0, '-', label="x_0", linewidth=1)
        plt.plot(ts, u_pc1, '-', label="x_1", linewidth=1)
        plt.plot(ts, u_pc2, '-', label="x_2", linewidth=1)
        plt.xlabel('time')
        plt.ylabel('displacement')
        plt.legend()
        plt.savefig('output/u_pc.png', format='png')
    # or save displacement
    u_pc_data=FileWriter('output/U_pc.out')
    for i in range(len(ts)) :
        u_pc_data.write("%f %f %f %f\n"%(ts[i],u_pc0[i],u_pc1[i],u_pc2[i]))
    u_pc_data.close()

