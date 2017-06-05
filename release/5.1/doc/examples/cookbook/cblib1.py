##############################################################################
#
# Copyright (c) 2009-2017 by The University of Queensland
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
from __future__ import division, print_function

__copyright__="""Copyright (c) 2009-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
A collection of routines to use in cookbook examples.
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

from esys.pycad import CurveLoop

# numpy for array handling
import numpy as np

# tools for dealing with PDEs - contains locator
from esys.escript.pdetools import Locator, Projector

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.weipa import saveVTK


# routine to find consecutive coordinates of a loop in pycad
def getLoopCoords(loop):
    # return all construction points of input
    temp = loop.getConstructionPoints()
    #create a numpy array for xyz components or construction points
    coords = np.zeros([len(temp),3],float)
    #place construction points in array
    for i in range(0,len(temp)):
        coords[i,:]=temp[i].getCoordinates()
    #return a numpy array
    return coords
    

##############################################################################
# subroutine: cbphones
# Allows us to record the values of a PDE at various 
# specified locations in the model.
# Arguments:
#   domain  : domain of model
#   U       : Current time state displacement solution.
#   phones  : Geophone Locations
#   dim     : model dimesions
#   savepath: where to output the data files local is default
##############################################################################
def cbphones(domain,U,phones,dim,savepath=""):
   #find the number of geophones
   nphones = len(phones)
   u_pot = np.zeros([nphones,dim],float)
   
   for i in range(0,nphones):
     # define the location of the phone source 
     L=Locator(domain,numpy.array(phones[i]))
     # find potential at point source.
     temp = L.getValue(U)
     for j in range(0,dim):
       u_pot[i,j]=temp[j]

   # open file to save displacement at point source
   return u_pot

##############################################################################
# subroutine: wavesolver2d
# Can solve a generic 2D wave propagation problem with a
# point source in a homogeneous medium.
# Arguments:
#   domain  : domain to solve over
#   h       : time step
#   tend    : end time
#   lam, mu : lames constants for domain
#   rho : density of domain
#   U0  : magnitude of source
#   xc  : source location in domain (Vector)
#   savepath: where to output the data files
##############################################################################
def wavesolver2d(domain,h,tend,lam,mu,rho,U0,xc,savepath,output="vtk"):
   from esys.escript.linearPDEs import LinearPDE
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
   print("src_radius = ",src_radius)

   dunit=numpy.array([0.,1.]) # defines direction of point source

   
   # ... set initial values ....
   n=0
   # initial value of displacement at point source is constant (U0=0.01)
   # for first two time steps
   u=U0*(cos(length(x-xc)*3.1415/src_radius)+1)*whereNegative(length(x-xc)-src_radius)*dunit
   u_m1=u
   t=0

   u_pot = cbphones(domain,u,[[0,500],[250,500],[400,500]],2)
   u_pc_x1 = u_pot[0,0]
   u_pc_y1 = u_pot[0,1]
   u_pc_x2 = u_pot[1,0]
   u_pc_y2 = u_pot[1,1]
   u_pc_x3 = u_pot[2,0]
   u_pc_y3 = u_pot[2,1]

   # open file to save displacement at point source
   u_pc_data=open(os.path.join(savepath,'U_pc.out'),'w')
   u_pc_data.write("%f %f %f %f %f %f %f\n"%(t,u_pc_x1,u_pc_y1,u_pc_x2,u_pc_y2,u_pc_x3,u_pc_y3))
 
#   while t<tend:
   while t<1.:
   
     # ... get current stress ....
      t=1.
##OLD WAY
      break
      g=grad(u)
      stress=lam*trace(g)*kmat+mu*(g+transpose(g))
     ### ... get new acceleration ....
     #mypde.setValue(X=-stress)          
     #a=mypde.getSolution()
     ### ... get new displacement ...
     #u_p1=2*u-u_m1+h*h*a
###NEW WAY
      mypde.setValue(X=-stress*(h*h),Y=(rho*2*u-rho*u_m1))
      u_p1 = mypde.getSolution()
     # ... shift displacements ....
      u_m1=u
      u=u_p1
     #stress = 
      t+=h
      n+=1
      print(n,"-th time step t ",t)
      u_pot = cbphones(domain,u,[[300.,200.],[500.,200.],[750.,200.]],2)

#     print "u at point charge=",u_pc
      u_pc_x1 = u_pot[0,0]
      u_pc_y1 = u_pot[0,1]
      u_pc_x2 = u_pot[1,0]
      u_pc_y2 = u_pot[1,1]
      u_pc_x3 = u_pot[2,0]
      u_pc_y3 = u_pot[2,1]
          
     # save displacements at point source to file for t > 0
      u_pc_data.write("%f %f %f %f %f %f %f\n"%(t,u_pc_x1,u_pc_y1,u_pc_x2,u_pc_y2,u_pc_x3,u_pc_y3))
 
     # ... save current acceleration in units of gravity and displacements 
     #saveVTK(os.path.join(savepath,"usoln.%i.vtu"%n),acceleration=length(a)/9.81,
     #displacement = length(u), tensor = stress, Ux = u[0] )
      if output == "vtk":
         saveVTK(os.path.join(savepath,"tonysol.%i.vtu"%n),output1 = length(u),tensor=stress)
      else:
         quT=qu.toListOfTuples()
    #Projector is used to smooth the data.
         proj=Projector(mymesh)
         smthT=proj(T)

    #move data to a regular grid for plotting
         xi,yi,zi = toRegGrid(smthT,mymesh,200,200,width,depth)

    # contour the gridded data, 
    # select colour
         pl.matplotlib.pyplot.autumn()
         pl.clf()
    # contour temperature
         CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
    # labels and formatting
         pl.clabel(CS, inline=1, fontsize=8)
         pl.title("Heat Refraction across a clinal structure.")
         pl.xlabel("Horizontal Displacement (m)")
         pl.ylabel("Depth (m)")
         pl.legend()
         if getMPIRankWorld() == 0: #check for MPI processing
            pl.savefig(os.path.join(saved_path,"heatrefraction001_cont.png"))

   u_pc_data.close()
   

##############################################################################
# subroutine: wavesolver2d
# Can solve a generic 2D wave propagation problem with a
# point source in a homogeneous medium with friction.
# Arguments:
#   domain  : domain to solve over
#   h       : time step
#   tend    : end time
#   lam, mu : lames constants for domain
#   rho : density of domain
#   U0  : magnitude of source
#   xc  : source location in domain (Vector)
#   savepath: where to output the data files
##############################################################################
def wavesolver2df(domain,h,tend,lam,mu,rho,U0,xc,savepath):
   x=domain.getX()
   # ... open new PDE ...
   mypde=LinearPDE(domain)
   #mypde.setSolverMethod(LinearPDE.LUMPING)
   mypde.setSymmetryOn()
   kmat = kronecker(domain)
   mypde.setValue(D=kmat)
   b=0.9

   # define small radius around point xc
   # Lsup(x) returns the maximum value of the argument x
   src_radius = 50#2*Lsup(domain.getSize())
   print("src_radius = ",src_radius)

   dunit=numpy.array([0.,1.]) # defines direction of point source

   
   # ... set initial values ....
   n=0
   # initial value of displacement at point source is constant (U0=0.01)
   # for first two time steps
   u=U0*(cos(length(x-xc)*3.1415/src_radius)+1)*whereNegative(length(x-xc)-src_radius)*dunit
   u_m1=u
   t=0

   u_pot = cbphones(domain,u,[[0,500],[250,500],[400,500]],2)
   u_pc_x1 = u_pot[0,0]
   u_pc_y1 = u_pot[0,1]
   u_pc_x2 = u_pot[1,0]
   u_pc_y2 = u_pot[1,1]
   u_pc_x3 = u_pot[2,0]
   u_pc_y3 = u_pot[2,1]

   # open file to save displacement at point source
   u_pc_data=open(os.path.join(savepath,'U_pc.out'),'w')
   u_pc_data.write("%f %f %f %f %f %f %f\n"%(t,u_pc_x1,u_pc_y1,u_pc_x2,u_pc_y2,u_pc_x3,u_pc_y3))
 
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
     y = ((rho/(-rho-b*h))*(u_m1-2*u))+(((b*h)/(-rho-(b*h)))*-u)
     mypde.setValue(X=-stress*((h*h)/(-rho-h*b)),Y=y)
     u_p1 = mypde.getSolution()
     # ... shift displacements ....
     u_m1=u
     u=u_p1
     #stress = 
     t+=h
     n+=1
     print(n,"-th time step t ",t)
     u_pot = cbphones(domain,u,[[300.,200.],[500.,200.],[750.,200.]],2)

#     print "u at point charge=",u_pc
     u_pc_x1 = u_pot[0,0]
     u_pc_y1 = u_pot[0,1]
     u_pc_x2 = u_pot[1,0]
     u_pc_y2 = u_pot[1,1]
     u_pc_x3 = u_pot[2,0]
     u_pc_y3 = u_pot[2,1]
           
     # save displacements at point source to file for t > 0
     u_pc_data.write("%f %f %f %f %f %f %f\n"%(t,u_pc_x1,u_pc_y1,u_pc_x2,u_pc_y2,u_pc_x3,u_pc_y3))
 
     # ... save current acceleration in units of gravity and displacements 
     #saveVTK(os.path.join(savepath,"usoln.%i.vtu"%n),acceleration=length(a)/9.81,
     #displacement = length(u), tensor = stress, Ux = u[0] )
     saveVTK(os.path.join(savepath,"tonysol.%i.vtu"%n),output1 = length(u),tensor=stress)

   u_pc_data.close()

