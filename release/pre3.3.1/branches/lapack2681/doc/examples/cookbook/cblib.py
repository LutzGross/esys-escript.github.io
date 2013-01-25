########################################################
#
# Copyright (c) 2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

from esys.pycad import CurveLoop

# numpy for array handling
import numpy as np
import pylab as pl
# tools for dealing with PDEs - contains locator
from esys.escript.pdetools import Locator, Projector

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE

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
	
# Calculate the location of quivers for a matplotlib plot
# quivshape :: [x,y] :: number of quivers in x and y direction
# lenxax :: length of model along x
# lenyax :: length of model along y
# qu :: gradient of escript solution ie grad(T)
def toQuivLocs(quivshape,lenxax,lenyax,qu):
    numquiv = quivshape[0]*quivshape[1] # total number of quivers
    dx = lenxax/quivshape[0]+1. # quiver x spacing
    dy = lenyax/quivshape[1]+1. # quiver y spacing
    qulocs=np.zeros([numquiv,2],float) # memory for quiver locations
    # fill qulocs
    for i in range(0,quivshape[0]-1):
    	for j in range(0,quivshape[1]-1):
			qulocs[i*quivshape[0]+j,:] = [dx+dx*i,dy+dy*j]
    # retreive values for quivers direction and shape from qu
    quL = Locator(qu.getFunctionSpace(),qulocs.tolist())
    qu = quL(qu) #array of dx,dy for quivers
    qu = np.array(qu) #turn into a numpy array
    qulocs = quL.getX() #returns actual locations from data
    qulocs = np.array(qulocs) #turn into a numpy array
    return qu,qulocs
    

 	

# Extract the X and Y coordinates of an array
# coords :: escript coordiantes from .getX function
def toXYTuple(coords):
    coords = np.array(coords.toListOfTuples()) #convert to Tuple
    coordX = coords[:,0]; coordY = coords[:,1] #X and Y components.
    return coordX,coordY

def toRegGrid(grid,domain,newx,newy,width,depth):
	oldspacecoords=domain.getX()
	gridT = grid.toListOfTuples(scalarastuple=False)
	cxy = Data(oldspacecoords, grid.getFunctionSpace())
	cX, cY = toXYTuple(cxy)
	xi = np.linspace(0.0,width,newx)	
	yi = np.linspace(depth,0.0,newy)
	# grid the data.
	zi = pl.matplotlib.mlab.griddata(cX,cY,gridT,xi,yi)
	return xi,yi,zi
	
def gradtoRegGrid(grid,domain,newx,newy,width,depth,dim):
	cxy = grid.getFunctionSpace().getX()
	gridT = grid.toListOfTuples()#(scalarastuple=False)
	#cxy = Data(oldspacecoords, grid.getFunctionSpace())
	cX, cY = toXYTuple(cxy)
	xi = np.linspace(0.0,width,newx)	
	yi = np.linspace(depth,0.0,newy)
	
	gridT = np.array(gridT)
	gridT = gridT[:,dim]
	
    # grid the data.
	
	zi = pl.matplotlib.mlab.griddata(cX,cY,gridT,xi,yi)
	return xi,yi,zi
########################################################
# subroutine: cbphones
# Allows us to record the values of a PDE at various 
# specified locations in the model.
# Arguments:
#	domain  : domain of model
#	U       : Current time state displacement solution.
#	phones  : Geophone Locations
#	dim     : model dimesions
#	savepath: where to output the data files local is default
########################################################
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
   print "src_radius = ",src_radius

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
     mypde.setValue(X=-stress*(h*h),Y=(rho*2*u-rho*u_m1))
     u_p1 = mypde.getSolution()
     # ... shift displacements ....
     u_m1=u
     u=u_p1
     #stress = 
     t+=h
     n+=1
     print n,"-th time step t ",t
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
   

########################################################
# subroutine: wavesolver2d
# Can solve a generic 2D wave propagation problem with a
# point source in a homogeneous medium with friction.
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
   print "src_radius = ",src_radius

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
     print n,"-th time step t ",t
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

# joel wrote this to create directories for you
import os
def needdirs(dirlist):
    for name in dirlist:
	if name == '':
		continue
	if not os.path.exists(name):
    	   try:
		os.makedirs(name)
    	   except OSError:
		if not os.path.exists(save_path):
	   	    raise
