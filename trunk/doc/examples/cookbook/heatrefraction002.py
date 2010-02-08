
########################################################
#
# Copyright (c) 2009-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# heatrefraction002.py
# Model steady state temperature distribution in three block model, 
# mesh from heatrefraction_mesher002.py 

#######################################################EXTERNAL MODULES
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
import matplotlib

#For interactive use, you can comment out the next line
matplotlib.use('agg') #It's just here for automated testing

# pylab for matplotlib and plotting
import pylab as pl
# cblib functions
from cblib import toXYTuple
from cblib2 import toQuivLocs

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

#################################################ESTABLISHING VARIABLES
qin=300.*Milli*W/(m*m) #our heat source temperature is now zero
Ti=290.15*K # Kelvin #the starting temperature of our iron bar
width=5000.0*m
depth=-6000.0*m

# the folder to gett our outputs from, leave blank "" for script path - 
# note these depen. are generated from heatrefraction_mesher001.py
saved_path = save_path= os.path.join("data","heatrefrac002")
mkDir(saved_path)

################################################ESTABLISHING PARAMETERS
## DOMAIN
mymesh=ReadMesh(os.path.join(saved_path,"heatrefraction_mesh003.fly"))
tpg = np.loadtxt(os.path.join(saved_path,"toppg"))
tpgx = tpg[:,0]
tpgy = tpg[:,1]
bpgl = np.loadtxt(os.path.join(saved_path,"botpgl"))
bpglx = bpgl[:,0]
bpgly = bpgl[:,1]
bpgr = np.loadtxt(os.path.join(saved_path,"botpgr"))
bpgrx = bpgr[:,0]
bpgry = bpgr[:,1]

# set up kappa (thermal conductivity across domain) using tags
kappa=Scalar(0,Function(mymesh))
kappa.setTaggedValue("top",2.0)
kappa.setTaggedValue("bottomleft",4.0)
kappa.setTaggedValue("bottomright",6.0)

#... generate functionspace...
#... open PDE ...
mypde=LinearPDE(mymesh)
mypde.setSymmetryOn()
#define first coefficient
mypde.setValue(A=kappa*kronecker(mymesh))

# ... set initial temperature ....
x=mymesh.getX()

qH=Scalar(0,FunctionOnBoundary(mymesh))
qH.setTaggedValue("linebottom",qin)
mypde.setValue(q=whereZero(x[1]),r=Ti)
mypde.setValue(y=qH)

###########################################################GET SOLUTION
T=mypde.getSolution()

##########################################################VISUALISATION
# calculate gradient of solution for quiver plot
qu=-kappa*grad(T)

# rearrage mymesh to suit solution function space      
oldspacecoords=mymesh.getX()
coords=Data(oldspacecoords, T.getFunctionSpace())

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
# contour the gridded data, plotting dots at the randomly 
# spaced data points.

pl.matplotlib.pyplot.autumn()              
CKL = pl.fill(tpgx,tpgy,'brown',label='2 W/m/k',zorder=-1000)
CKM = pl.fill(bpglx,bpgly,'red',label='4 W/m/k',zorder=-1000)
CKN = pl.fill(bpgrx,bpgry,'orange',label='6 W/m/k',zorder=-1000)              
              
CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
pl.clabel(CS, inline=1, fontsize=8)
pl.title("Heat Refraction across an anisotropic structure.")
pl.xlabel("Horizontal Displacement (m)")
pl.ylabel("Depth (m)")
pl.legend()
if getMPIRankWorld() == 0: #check for MPI processing
	pl.savefig(os.path.join(saved_path,"heatrefraction002_cont.png"))

QUIV=pl.quiver(qulocs[:,0],qulocs[:,1],qu[:,0],qu[:,1],\
                angles='xy',color="white")
pl.title("Heat Refraction across an anisotropic structure\
            \n with gradient quivers.")
if getMPIRankWorld() == 0: #check for MPI processing
	pl.savefig(os.path.join(saved_path,"heatrefraction002_contqu.png"))
