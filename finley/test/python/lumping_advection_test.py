

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

   a simple comparison for row-sum and HRZ lumping in case of the advection equation

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



v=numpy.array([-1,0])
n=5.

def ref_u(x,t):
   return whereNonPositive(x[0]-t)
   #return exp(-5*(x[0]-t)**2)

def runTaylorGalerkinIncremental(order):
   domain=finley.Rectangle(100,10,order)
   x=domain.getX()
   
   
   # test Velet scheme 
   dt=inf(domain.getSize()/length(v))*(1./6.)
   q=whereZero(x[0])+whereZero(x[0]-1.)
   
   mypde_f=LinearSinglePDE(domain)
   mypde_f.setSymmetryOn()
   mypde_f.setValue(D=1,q=q)
   u_f=ref_u(x,0)
   
   mypde_HRZ=LinearSinglePDE(domain)
   mypde_HRZ.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
   mypde_HRZ.setValue(D=1,q=q)
   u_HRZ=ref_u(x,0)
   
   mypde_RS=LinearSinglePDE(domain)
   mypde_RS.getSolverOptions().setSolverMethod(SolverOptions.ROWSUM_LUMPING)
   mypde_RS.setValue(D=1,q=q)
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
   
   while t< 1./Lsup(v):
       t+=dt
       u=ref_u(x,t)
       mypde_f.setValue(X=-dt/2.*u_f*v, r=ref_u(x,t-dt/2)-u_f)
       mypde_HRZ.setValue(X=-dt/2.*u_HRZ*v, r=ref_u(x,t-dt/2)-u_f)
       mypde_RS.setValue(X=-dt/2.*u_RS*v, r=ref_u(x,t-dt/2)-u_f)
    
       u_f_h=u_f+mypde_f.getSolution()
       u_HRZ_h=u_HRZ+mypde_HRZ.getSolution()
       u_RS_h=u_RS+mypde_RS.getSolution()
   
       mypde_f.setValue(X=-dt*u_f_h*v, r=u-u_f)
       mypde_HRZ.setValue(X=-dt*u_HRZ_h*v, r=u-u_HRZ)
       mypde_RS.setValue(X=-dt*u_RS_h*v, r=u-u_RS)
    
       u_f=u_f+mypde_f.getSolution()
       u_HRZ=u_HRZ+mypde_HRZ.getSolution()
       u_RS=u_RS+mypde_RS.getSolution()

       t_list.append(t)
       u_list.append(l(u))
       f_list.append(l(u_f))
       HRZ_list.append(l(u_HRZ))
       RS_list.append(l(u_RS))
       print(t_list[-1], u_list[-1], f_list[-1], HRZ_list[-1] , RS_list[-1], " : ",sup(u))
     
   
   import matplotlib.pyplot as plt
   if getMPIRankWorld() == 0:
         plt.clf()
         plt.plot(t_list, u_list, '-', label="exact", linewidth=1)
         plt.plot(t_list, f_list, '-', label="full", linewidth=1)
         plt.plot(t_list, HRZ_list, '-', label="HRZ lumping", linewidth=1)
         plt.plot(t_list, RS_list, '-', label="row sum lumping", linewidth=1)
         plt.axis([0.,max(t_list),-.3,2.])
         plt.xlabel('time')
         plt.ylabel('displacement')
         plt.legend()
         plt.savefig('lumping_SUPG_du_%d.png'%order, format='png')

def runTaylorGalerkinDirect(order):
   domain=finley.Rectangle(500,10,order)
   x=domain.getX()
   
   
   # test Velet scheme 
   dt=inf(domain.getSize()/length(v))*(1./6.)
   q=whereZero(x[0])+whereZero(x[0]-1.)
   
   mypde_f=LinearSinglePDE(domain)
   mypde_f.setSymmetryOn()
   mypde_f.setValue(D=1,q=q)
   u_f=ref_u(x,0)
   
   mypde_HRZ=LinearSinglePDE(domain)
   mypde_HRZ.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
   mypde_HRZ.setValue(D=1,q=q)
   u_HRZ=ref_u(x,0)
   
   mypde_RS=LinearSinglePDE(domain)
   mypde_RS.getSolverOptions().setSolverMethod(SolverOptions.ROWSUM_LUMPING)
   mypde_RS.setValue(D=1,q=q)
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
   
   while t< 1./Lsup(v):
       t+=dt
       u=ref_u(x,t)
       mypde_f.setValue(X=-dt/2.*u_f*v, r=ref_u(x,t-dt/2), Y=u_f )
       mypde_HRZ.setValue(X=-dt/2.*u_HRZ*v, r=ref_u(x,t-dt/2), Y=u_HRZ)
       mypde_RS.setValue(X=-dt/2.*u_RS*v, r=ref_u(x,t-dt/2), Y= u_RS)
    
       u_f_h=mypde_f.getSolution()
       u_HRZ_h=mypde_HRZ.getSolution()
       u_RS_h=mypde_RS.getSolution()
   
       mypde_f.setValue(X=-dt*u_f_h*v, r=u, Y=u_f )
       mypde_HRZ.setValue(X=-dt*u_HRZ_h*v, r=u, Y=u_HRZ)
       mypde_RS.setValue(X=-dt*u_RS_h*v, r=u, Y= u_RS)
    
       u_f=mypde_f.getSolution()
       u_HRZ=mypde_HRZ.getSolution()
       u_RS=mypde_RS.getSolution()

       t_list.append(t)
       u_list.append(l(u))
       f_list.append(l(u_f))
       HRZ_list.append(l(u_HRZ))
       RS_list.append(l(u_RS))
       print(t_list[-1], u_list[-1], f_list[-1], HRZ_list[-1] , RS_list[-1], " : ",sup(u))
     
   
   import matplotlib.pyplot as plt
   if getMPIRankWorld() == 0:
         plt.clf()
         plt.plot(t_list, u_list, '-', label="exact", linewidth=1)
         plt.plot(t_list, f_list, '-', label="full", linewidth=1)
         plt.plot(t_list, HRZ_list, '-', label="HRZ lumping", linewidth=1)
         plt.plot(t_list, RS_list, '-', label="row sum lumping", linewidth=1)
         plt.axis([0.,max(t_list),-.3,2.])
         plt.xlabel('time')
         plt.ylabel('displacement')
         plt.legend()
         plt.savefig('lumping_SUPG_u_%d.png'%order, format='png')


# runTaylorGalerkinIncremental(1)
# runTaylorGalerkinIncremental(2)
runTaylorGalerkinDirect(1)
# runTaylorGalerkinDirect(2)

