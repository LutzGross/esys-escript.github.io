
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

############################################################FILE HEADER
# twodheatdiff002.py
# Model temperature diffusion between a granite intrusion and sandstone 
# country rock. This is a two dimensional problem with the granite as a
# heat source.

#######################################################EXTERNAL MODULES
#To solve the problem it is necessary to import the modules we require.
#This imports everything from the escript library
from esys.escript import *
# This defines the LinearPDE module as LinearPDE
from esys.escript.linearPDEs import LinearPDE 
# This imports the rectangle domain function from finley.
from esys.finley import Rectangle 
# A useful unit handling package which will make sure all our units
# match up in the equations under SI.
from esys.escript.unitsSI import *
#For interactive use, you can comment out the next two lines
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
import pylab as pl #Plotting package.
import numpy as np #Array package.
import os #This package is necessary to handle saving our data.
from cblib import toXYTuple, needdirs

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

#################################################ESTABLISHING VARIABLES
#PDE related
mx = 600*m #meters - model length
my = 600*m #meters - model width
ndx = 100 #mesh steps in x direction 
ndy = 100 #mesh steps in y direction
r = 200*m #meters - radius of intrusion
ic = [300, 0] #centre of intrusion (meters)
q=0.*Celsius #our heat source temperature is now zero

## Intrusion Variables - Granite
Ti=2273.*Celsius # Kelvin -the starting temperature of our RHS Block
rhoi = 2750*kg/m**3 #kg/m^{3} density of granite
cpi = 790.*J/(kg*K) #j/Kg.K thermal capacity
rhocpi = rhoi*cpi	#DENSITY * SPECIFIC HEAT
eta=0.  # RADIATION CONDITION
kappai=2.2*W/m/K #watts/m.K thermal conductivity
## Country Rock Variables - Sandstone
Tc = 473*Celsius # Kelvin #the starting temperature of our country rock
rhoc = 2000*kg/m**3 #kg/m^{3} density
cpc = 920.*J/(kg*K) #j/kg.k specific heat
rhocpc = rhoc*cpc #DENSITY * SPECIFIC HEAT
kappac = 1.9*W/m/K #watts/m.K thermal conductivity

#Script/Iteration Related
t=0. #our start time, usually zero
tday=100*365. #the time we want to end the simulation in days
tend=tday*24*60*60
outputs = 200 # number of time steps required.
h=(tend-t)/outputs #size of time step
#user warning
print "Expected Number of Output Files is: ", outputs
print "Step size is: ", h/(24.*60*60), "days"
i=0 #loop counter 
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","twodheatdiff")
needdirs([save_path])
########## note this folder path must exist to work ###################

################################################ESTABLISHING PARAMETERS
#generate domain using rectangle
model = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
#extract finite points - the solution points
x=model.getX()
#create the PDE
mypde=LinearPDE(model) #assigns a domain to our PDE
mypde.setSymmetryOn() #set the fast solver on for symmetry
#establish location of boundary between two materials
bound = length(x-ic)-r #where the boundary will be located
A = (kappai)*whereNegative(bound)+(kappac)*wherePositive(bound)
D = (rhocpi/h)*whereNegative(bound)+(rhocpc/h)*wherePositive(bound)
#define our PDE coeffs
mypde.setValue(A=A*kronecker(model),D=D,d=eta,y=eta*Tc)
#set initial temperature
T= Ti*whereNegative(bound)+Tc*wherePositive(bound)

# rearrage mymesh to suit solution function space for contouring      
oldspacecoords=model.getX()
coords=Data(oldspacecoords, T.getFunctionSpace())
#coords = np.array(coords.toListOfTuples())
coordX, coordY = toXYTuple(coords)
# create regular grid
xi = np.linspace(0.0,mx,100)
yi = np.linspace(0.0,my,100)

########################################################START ITERATION
while t<=tend:
      i+=1 #counter
      t+=h #curretn time
      Y = T*D #
      mypde.setValue(Y=Y)
      T=mypde.getSolution()
      tempT = T.toListOfTuples(scalarastuple=False)
      # grid the data.
      zi = pl.matplotlib.mlab.griddata(coordX,coordY,tempT,xi,yi)
      # contour the gridded data, plotting dots at the 
      # randomly spaced data points.
      pl.matplotlib.pyplot.autumn()
      pl.contourf(xi,yi,zi,10)
      CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
      pl.clabel(CS, inline=1, fontsize=8)
      pl.axis([0,600,0,600])
      pl.title("Heat diffusion from an intrusion.")
      pl.xlabel("Horizontal Displacement (m)")
      pl.ylabel("Depth (m)")
      if getMPIRankWorld() == 0:
	      pl.savefig(os.path.join(save_path,\
	                               "heatrefraction%03d.png"%i))
      pl.clf()            

# compile the *.png files to create an *.avi video that shows T change
# with time. This opperation uses linux mencoder.
os.system("mencoder mf://"+save_path+"/*.png -mf type=png:\
w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
twodheatdiff001tempT.avi")
