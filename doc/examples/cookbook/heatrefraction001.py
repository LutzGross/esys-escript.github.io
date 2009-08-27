
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

# To solve the problem it is necessary to import the modules we 
# require.
# This imports everything from the escript library
from esys.escript import * 
# This defines LinearPDE as LinearPDE
from esys.escript.linearPDEs import LinearPDE, Poisson 
# This imports the rectangle domain function from finley
from esys.finley import Rectangle, ReadMesh, Domain 
# This package is necessary to handle saving our data.
import os
# A useful unit handling package which will make sure all our units
# match up in the equations.
from esys.escript.unitsSI import * 
# numpy for array handling
import numpy as np
# pylab for matplotlib and plotting
import pylab as pl
# cblib functions
from cblib import toQuivLocs, toXYTuple

##ESTABLISHING VARIABLES
qin=70*Milli*W/(m*m) #our heat source temperature is now zero
Ti=290.15*K # Kelvin #the starting temperature of our iron bar
width=5000.0*m
depth=-6000.0*m

# the folder to gett our outputs from, leave blank "" for script path - 
# note these depen. are generated from heatrefraction_mesher001.py
saved_path = "data/heatrefrac001" 

###### 2 BLOCK MODEL #########
## DOMAIN
## Anticline
mymesh=ReadMesh(os.path.join(saved_path,"heatrefraction_mesh001.fly"))
tpg = np.loadtxt(os.path.join(saved_path,"toppg"))
tpgx = tpg[:,0]
tpgy = tpg[:,1]
bpg = np.loadtxt(os.path.join(saved_path,"botpg"))
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

qH=Scalar(0,FunctionOnBoundary(mymesh))
qH.setTaggedValue("linebottom",qin)
mypde.setValue(q=whereZero(x[1]),r=Ti)
mypde.setValue(y=qH)#,r=17*Celsius)

# get steady state solution and export to vtk.
T=mypde.getSolution()
#saveVTK("tempheatrefract.xml",sol=T, q=-kappa*grad(T))

# rearrage mymesh to suit solution function space      
oldspacecoords=mymesh.getX()
coords=Data(oldspacecoords, T.getFunctionSpace())

# calculate gradient of solution for quiver plot
qu=-kappa*grad(T)
quivshape = [20,20] #quivers x and quivers y
# function to calculate quiver locations
qu,qulocs = toQuivLocs(quivshape,width,depth,qu)

kappaT = kappa.toListOfTuples(scalarastuple=False)
coordsK = Data(oldspacecoords, kappa.getFunctionSpace())
coordKX, coordKY = toXYTuple(coordsK)
      
tempT = T.toListOfTuples(scalarastuple=False)
coordX, coordY = toXYTuple(coords)

xi = np.linspace(0.0,width,100)
yi = np.linspace(depth,0.0,100)
# grid the data.
zi = pl.matplotlib.mlab.griddata(coordX,coordY,tempT,xi,yi)
ziK = pl.matplotlib.mlab.griddata(coordKX,coordKY,kappaT,xi,yi)
# contour the gridded data, plotting dots at the randomly spaced data points.

pl.matplotlib.pyplot.autumn()
CKL = pl.fill(tpgx,tpgy,'brown',bpgx,bpgy,'red',zorder=-1000)
#~ CK = pl.contourf(xi,yi,ziK,2)
CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
pl.clabel(CS, inline=1, fontsize=8)
pl.title("Heat Refraction across a clinal structure.")
pl.xlabel("Horizontal Displacement (m)")
pl.ylabel("Depth (m)")
#~ CB = pl.colorbar(CS, shrink=0.8, extend='both')
pl.savefig(os.path.join(saved_path,"heatrefraction001_cont.png"))

QUIV=pl.quiver(qulocs[:,0],qulocs[:,1],qu[:,0],qu[:,1],angles='xy',color="white")
pl.title("Heat Refraction across a clinal structure \n with gradient quivers.")
pl.savefig(os.path.join(saved_path,"heatrefraction001_contqu.png"))
