
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

"""
   a simple comparison for row-sum and HRZ lumping in case of the wave equation
"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator
from esys import finley
from math import pi



c=1.
n=5.

def ref_u(x,t):
   return sin(pi*n*(x[0]-c*t))

def runValetAcceleration(order):
   domain=finley.Rectangle(100,10,order)
   x=domain.getX()
   
   
   # test Velet scheme 
   dt=inf(domain.getSize()/c)*(1./6.)
   q=whereZero(x[0])+whereZero(x[0]-1.)
   
   mypde_f=LinearSinglePDE(domain)
   mypde_f.setSymmetryOn()
   mypde_f.setValue(D=1,q=q)
   u_f_old=ref_u(x,-dt)
   u_f=ref_u(x,0)
   
   mypde_HRZ=LinearSinglePDE(domain)
   mypde_HRZ.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
   mypde_HRZ.setValue(D=1,q=q)
   u_HRZ_old=ref_u(x,-dt)
   u_HRZ=ref_u(x,0)
   
   mypde_RS=LinearSinglePDE(domain)
   mypde_RS.getSolverOptions().setSolverMethod(SolverOptions.ROWSUM_LUMPING)
   mypde_RS.setValue(D=1,q=q)
   u_RS_old=ref_u(x,-dt)
   u_RS=ref_u(x,0)
   
   l=Locator(domain,[0.5,0.5])
   t=0
   
   u=ref_u(x,t)
   t_list=[t]
   u_list=[l(u)]
   f_list=[l(u_f)]
   HRZ_list=[l(u_HRZ)]
   RS_list=[l(u_RS)]
   print(t_list[-1], u_list[-1], f_list[-1], HRZ_list[-1] , RS_list[-1])
   
   while t< 4/n/c:
       t+=dt
       u=ref_u(x,t)
       mypde_f.setValue(X=-c**2*grad(u_f), r=-c**2*u)
       mypde_HRZ.setValue(X=-c**2*grad(u_HRZ), r=-c**2*u)
       mypde_RS.setValue(X=-c**2*grad(u_RS), r=-c**2*u)
    
       a_f=mypde_f.getSolution()
       a_HRZ=mypde_HRZ.getSolution()
       a_RS=mypde_RS.getSolution()
   
       u_f, u_f_old = 2*u_f-u_f_old + dt**2*a_f ,  u_f
       u_HRZ, u_HRZ_old = 2*u_HRZ-u_HRZ_old + dt**2*a_HRZ ,  u_HRZ
       u_RS, u_RS_old = 2*u_RS-u_RS_old + dt**2*a_RS ,  u_RS
      
       t_list.append(t)
       u_list.append(l(u))
       f_list.append(l(u_f))
       HRZ_list.append(l(u_HRZ))
       RS_list.append(l(u_RS))
       print(t_list[-1], u_list[-1], f_list[-1], HRZ_list[-1] , RS_list[-1])
     
   
   import matplotlib.pyplot as plt
   if getMPIRankWorld() == 0:
         plt.clf()
         plt.plot(t_list, u_list, '-', label="exact", linewidth=1)
         plt.plot(t_list, f_list, '-', label="full", linewidth=1)
         plt.plot(t_list, HRZ_list, '-', label="HRZ lumping", linewidth=1)
         plt.plot(t_list, RS_list, '-', label="row sum lumping", linewidth=1)
         plt.axis([0.,max(t_list),-1.3,1.3])
         plt.xlabel('time')
         plt.ylabel('displacement')
         plt.legend()
         plt.savefig('lumping_valet_a_%d.png'%order, format='png')

def runValetDisplacement(order):
   domain=finley.Rectangle(100,10,order)
   x=domain.getX()
   
   
   # test Velet scheme 
   dt=inf(domain.getSize()/c)*(1./6.)
   q=whereZero(x[0])+whereZero(x[0]-1.)
   
   mypde_f=LinearSinglePDE(domain)
   mypde_f.setSymmetryOn()
   mypde_f.setValue(D=1,q=q)
   u_f_old=ref_u(x,-dt)
   u_f=ref_u(x,0)
   
   mypde_HRZ=LinearSinglePDE(domain)
   mypde_HRZ.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
   mypde_HRZ.setValue(D=1,q=q)
   u_HRZ_old=ref_u(x,-dt)
   u_HRZ=ref_u(x,0)
   
   mypde_RS=LinearSinglePDE(domain)
   mypde_RS.getSolverOptions().setSolverMethod(SolverOptions.ROWSUM_LUMPING)
   mypde_RS.setValue(D=1,q=q)
   u_RS_old=ref_u(x,-dt)
   u_RS=ref_u(x,0)
   
   l=Locator(domain,[0.5,0.5])
   t=0
   
   u=ref_u(x,t)
   t_list=[t]
   u_list=[l(u)]
   f_list=[l(u_f)]
   HRZ_list=[l(u_HRZ)]
   RS_list=[l(u_RS)]
   print(t_list[-1], u_list[-1], f_list[-1], HRZ_list[-1] , RS_list[-1])
   
   while t< 4/n/c:
       t+=dt
       u=ref_u(x,t)
       mypde_f.setValue(X=-(c*dt)**2*grad(u_f), r=u, Y=2*u_f-u_f_old )
       mypde_HRZ.setValue(X=-(c*dt)**2*grad(u_HRZ), r=u, Y= 2*u_HRZ-u_HRZ_old )
       mypde_RS.setValue(X=-(c*dt)**2*grad(u_RS), r=u, Y= 2*u_RS-u_RS_old )
    
       u_f, u_f_old =mypde_f.getSolution(), u_f
       u_HRZ, u_HRZ_old =mypde_HRZ.getSolution(), u_HRZ
       u_RS, u_RS_old =mypde_RS.getSolution(), u_RS
   
       t_list.append(t)
       u_list.append(l(u))
       f_list.append(l(u_f))
       HRZ_list.append(l(u_HRZ))
       RS_list.append(l(u_RS))
       print(t_list[-1], u_list[-1], f_list[-1], HRZ_list[-1] , RS_list[-1])
     
   
   import matplotlib.pyplot as plt
   if getMPIRankWorld() == 0:
         plt.clf()
         plt.plot(t_list, u_list, '-', label="exact", linewidth=1)
         plt.plot(t_list, f_list, '-', label="full", linewidth=1)
         plt.plot(t_list, HRZ_list, '-', label="HRZ lumping", linewidth=1)
         plt.plot(t_list, RS_list, '-', label="row sum lumping", linewidth=1)
         plt.axis([0.,max(t_list),-1.3,1.3])
         plt.xlabel('time')
         plt.ylabel('displacement')
         plt.legend()
         plt.savefig('lumping_valet_u_%d.png'%order, format='png')


#runValetAcceleration(1)
#runValetAcceleration(2)
runValetDisplacement(1)
runValetDisplacement(2)

