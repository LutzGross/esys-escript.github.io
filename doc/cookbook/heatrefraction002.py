
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.linearPDEs import LinearPDE, Poisson # This defines LinearPDE as LinearPDE
from esys.finley import Rectangle, ReadMesh, Domain # This imports the rectangle domain function from finley
import os #This package is necessary to handle saving our data.
from esys.escript.unitsSI import * # A useful unit handling package which will make sure all our units match up in the equations.
import numpy as np
import pylab as pl

from esys.escript.pdetools import *

##ESTABLISHING VARIABLES
#DOMAIN
mymesh = ReadMesh("heatrefraction_mesh003.fly")
#~ print Function(mymesh).getListOfTags()

qin=70*Milli*W/(m*m) #our heat source temperature is now zero
Ti=290.15*K # Kelvin #the starting temperature of our iron bar

# set up kappa (thermal conductivity across domain) using tags
kappa=Scalar(0,Function(mymesh))
kappa.setTaggedValue("tblockloop",2.0)
kappa.setTaggedValue("bblockloopl",3.0)
kappa.setTaggedValue("bblockloopr",4.0)

#... generate functionspace...
#... open PDE ...
mypde=LinearPDE(mymesh)
#define first coefficient
mypde.setValue(A=kappa*kronecker(mymesh))

# ... set initial temperature ....
x=mymesh.getX()

qH=qin*whereZero(x[1]+6000)
mypde.setValue(q=Ti*whereZero(x[1]),r=Ti)
mypde.setValue(Y=qH,y=17*Celsius)

# get steady state solution and export to vtk.
T=mypde.getSolution()
saveVTK("tempheatrefract.xml",sol=T, q=-kappa*grad(T))

# rearrage mymesh to suit solution function space      
oldspacecoords=mymesh.getX()
coords=Data(oldspacecoords, T.getFunctionSpace())

# calculate gradient of solution for quiver plot
qu=-kappa*grad(T)

# create quiver locations
quivshape = [20,20] #quivers x and quivers y
numquiv = quivshape[0]*quivshape[1] # total number of quivers
maxx = 5000. # maximum x displacement
dx = maxx/quivshape[0]+1. # quiver x spacing
maxy = -6000. # maximum y displacement
dy = maxy/quivshape[1]+1. # quiver y spacing
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

kappaT = kappa.toListOfTuples(scalarastuple=False)
coordsK = Data(oldspacecoords, kappa.getFunctionSpace())
coordsK = np.array(coordsK.toListOfTuples())
coordKX = coordsK[:,0]
coordKY = coordsK[:,1]
      
coords = np.array(coords.toListOfTuples())
tempT = T.toListOfTuples(scalarastuple=False)

coordX = coords[:,0]
coordY = coords[:,1]

xi = np.linspace(0.0,5000.0,100)
yi = np.linspace(-6000.0,0.0,100)
# grid the data.
zi = pl.matplotlib.mlab.griddata(coordX,coordY,tempT,xi,yi)
ziK = pl.matplotlib.mlab.griddata(coordKX,coordKY,kappaT,xi,yi)
# contour the gridded data, plotting dots at the randomly spaced data points.

pl.matplotlib.pyplot.autumn()
CKL = pl.contour(xi,yi,ziK,1,linewidths=1.0)
CK = pl.contourf(xi,yi,ziK,1)
CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
pl.clabel(CS, inline=1, fontsize=8)
CB = pl.colorbar(CS, shrink=0.8, extend='both')
pl.savefig("temp.png")

#~ dx = np.zeros(82,float)
#~ dy = np.zeros(82,float)
#~ vlocs = np.zeros([82,2],float)
#~ delx = 5000.0/100
#~ dely = 6000.0/100
#~ ind = 0
#~ 
#~ for i in range(1,10):
	#~ for j in range(1,10):
		#~ iloc = i*10
		#~ jloc = j*10
		#~ ind += 1
		#~ 
		#~ vlocs[ind,0] = iloc*delx; vlocs[ind,1]=-jloc*dely
		#~ 
		#~ print zi[iloc,jloc], zi[iloc+5,jloc]
		#~ print zi[iloc,jloc], zi[iloc,jloc+5]
		#~ dx[ind] = (zi[iloc,jloc] - zi[iloc+5,jloc])/delx
		#~ dy[ind] = (zi[iloc,jloc] - zi[iloc,jloc+5])/dely

#~ zpoints = [zi[iloc,jloc],zi[iloc+2,jloc],zi[iloc+4,jloc]]
#~ xpoints = [vlocs[ind,0],vlocs[ind+2,0],vlocs[ind+4,0]]
#~ print pl.matplotlib.mlab.slopes(xpoints,zpoints)
#~ dx[ind] = pl.matplotlib.mlab.slopes(xpoints,zpoints)
		#~ 
#~ ypoints = [vlocs[ind,1],vlocs[ind+2,1],vlocs[ind+4,]]
#~ zpoints = [zi[iloc,jloc],zi[iloc,jloc+2],zi[iloc,jloc+4]]
#~ dy[ind] = pl.matplotlib.mlab.slopes(ypoints,zpoints)
#~ 
pl.quiver(qulocs[:,0],qulocs[:,1],qu[:,0],qu[:,1],angles='xy',color="green")
pl.savefig("temp2.png")
