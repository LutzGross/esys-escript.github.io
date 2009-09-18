
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

############################################################FILE HEADER
# heatrefraction001.py
# Model steady state temperature distribution in two block model, mesh
# from heatrefraction_mesher001.py 

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we 
# require.



# This imports everything from the escript library
from esys.escript import * 
# This defines LinearPDE as LinearPDE
from esys.escript.linearPDEs import LinearPDE, Poisson 
from esys.escript.pdetools import Projector
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
from cblib import toQuivLocs, toXYTuple, needdirs, toRegGrid, gradtoRegGrid

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

#################################################ESTABLISHING VARIABLES
qin=70*Milli*W/(m*m) #our heat source temperature is now zero
Ti=290.15*K # Kelvin #the starting temperature of our iron bar
width=5000.0*m
depth=-6000.0*m

# the folder to gett our outputs from, leave blank "" for script path - 
# note these depen. are generated from heatrefraction_mesher001.py
saved_path = save_path= os.path.join("data","heatrefrac001" )
needdirs([saved_path])

################################################ESTABLISHING PARAMETERS
## DOMAIN
mymesh=ReadMesh(os.path.join(saved_path,"heatrefraction_mesh001.fly"))
tpg = np.loadtxt(os.path.join(saved_path,"toppg"))
tpgx = tpg[:,0]
tpgy = tpg[:,1]
bpg = np.loadtxt(os.path.join(saved_path,"botpg"))
bpgx = bpg[:,0]
bpgy = bpg[:,1]

# set up kappa (thermal conductivity across domain) using tags
kappa=Scalar(0,Function(mymesh))
kappa.setTaggedValue("top",2.0)
kappa.setTaggedValue("bottom",4.0)

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
mypde.setValue(y=qH)

###########################################################GET SOLUTION
T=mypde.getSolution()

##########################################################VISUALISATION
# calculate gradient of solution for quiver plot
qu=-kappa*grad(T)
quT=qu.toListOfTuples()

#Projector is used to smooth the data.
proj=Projector(mymesh)
smthT=proj(T)

quivshape = [20,20] #quivers x and quivers y
# function to calculate quiver locations
qu,qulocs = toQuivLocs(quivshape,width,depth,qu)

#move data to a regular grid for plotting
xi,yi,zi = toRegGrid(smthT,mymesh,200,200,width,depth)

# contour the gridded data, 
# plotting dots at the randomly spaced data points.

# select colour
pl.matplotlib.pyplot.autumn()
# plot polygons for boundaries
CKL = pl.fill(tpgx,tpgy,'brown',label='2 W/m/k',zorder=-1000)
CKM = pl.fill(bpgx,bpgy,'red',label='4 W/m/k',zorder=-1000)
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

#Quiver Plot qulocs -> tail location, qu -> quiver length/direction
QUIV=pl.quiver(qulocs[:,0],qulocs[:,1],qu[:,0],qu[:,1],\
                angles='xy',color="white")
pl.title("Heat Refraction across a clinal structure\n with\
gradient quivers.")
if getMPIRankWorld() == 0: #check for MPI processing
	pl.savefig(os.path.join(saved_path,"heatrefraction001_contqu.png"))

#Temperature Depth Profile along x[50]
cut=int(len(xi)/2)
pl.clf()
pl.plot(zi[:,cut],yi)
pl.title("Heat Refraction Temperature Depth Profile")
pl.xlabel("Temperature (K)")
pl.ylabel("Depth (m)")
if getMPIRankWorld() == 0: #check for MPI processing
    pl.savefig(os.path.join(saved_path,"heatrefraction001_tdp.png"))
    
#Temperature Gradient Profile along x[50]



pl.clf()
# grid the data.
qu=proj(-kappa*grad(T))
xiq,yiq,ziq = gradtoRegGrid(qu,mymesh,200,200,width,depth,1)
cut=int(len(xiq)/2)
pl.plot(ziq[:,cut]*1000.,yiq)
pl.title("Heat Flow Depth Profile")
pl.xlabel("Heat Flow (mW/m^2)")
pl.ylabel("Depth (m)")
if getMPIRankWorld() == 0: #check for MPI processing
	pl.savefig(os.path.join(saved_path,"heatrefraction001_hf.png"))

pl.clf()
zT=proj(-grad(T))

xt,yt,zt=gradtoRegGrid(zT,mymesh,200,200,width,depth,1)
cut=int(len(xt)/2)
pl.plot(zt[:,cut]*1000.,yt)
pl.title("Heat Refraction Temperature Gradient \n Depth Profile")
pl.xlabel("Temperature (K/Km)")
pl.ylabel("Depth (m)")
if getMPIRankWorld() == 0: #check for MPI processing
    pl.savefig(os.path.join(saved_path,"heatrefraction001_tgdp.png"))
    
pl.clf()
xk,yk,zk = toRegGrid(kappa,mymesh,200,200,width,depth)
cut=int(len(xk)/2)
pl.plot(zk[:,cut],yk)
pl.title("Heat Refraction Thermal Conductivity Depth Profile")
pl.xlabel("Conductivity (W/K/m)")
pl.ylabel("Depth (m)")
pl.axis([1,5,-6000,0])
if getMPIRankWorld() == 0: #check for MPI processing
    pl.savefig(os.path.join(saved_path,"heatrefraction001_tcdp.png"))
    
