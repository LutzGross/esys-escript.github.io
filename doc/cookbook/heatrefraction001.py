
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
qin=70*Milli*W/(m*m) #our heat source temperature is now zero
Ti=290.15*K # Kelvin #the starting temperature of our iron bar

###### 2 BLOCK MODEL #########
## DOMAIN
## Anticline
mymesh = ReadMesh("heatrefraction_mesh001.fly")
tpg = np.loadtxt("toppg")
tpgx = tpg[:,0]
tpgy = tpg[:,1]
bpg = np.loadtxt("botpg")
bpgx = bpg[:,0]
bpgy = bpg[:,1]
## Syncline
#~ mymesh = ReadMesh("heatrefraction_mesh002.fly")

# set up kappa (thermal conductivity across domain) using tags
kappa=Scalar(0,Function(mymesh))
kappa.setTaggedValue("top",2.0)
kappa.setTaggedValue("bottom",4.0)

###### 3 BLOCK MODEL #########
# set up kappa (thermal conductivity across domain) using tags
#~ mymesh = ReadMesh("heatrefraction_mesh003.fly")
#~ kappa=Scalar(0,Function(mymesh))
#~ kappa.setTaggedValue("top",2.0)
#~ kappa.setTaggedValue("bottomleft",3.0)
#~ kappa.setTaggedValue("bottomright",4.0)

#... generate functionspace...
#... open PDE ...
mypde=LinearPDE(mymesh)
#define first coefficient
mypde.setValue(A=kappa*kronecker(mymesh))

# ... set initial temperature ....
x=mymesh.getX()

qH=qin*whereZero(x[1]-inf(x[1]))
mypde.setValue(q=whereZero(x[1]),r=Ti)
mypde.setValue(y=qH)#,r=17*Celsius)

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
CKL = pl.fill(tpgx,tpgy,'brown',bpgx,bpgy,'red',zorder=-1000)
#~ CK = pl.contourf(xi,yi,ziK,2)
CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
pl.clabel(CS, inline=1, fontsize=8)
#~ CB = pl.colorbar(CS, shrink=0.8, extend='both')
pl.savefig("temp.png")

QUIV = pl.quiver(qulocs[:,0],qulocs[:,1],qu[:,0],qu[:,1],angles='xy',color="white")
pl.savefig("temp2.png")
