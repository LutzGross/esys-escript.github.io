
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
# onedheatdiff002.py
# Model temperature diffusion between two granite blocks. This is a one
# dimensional problem with no heat source and a single heat disparity.

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
import pylab as pl #Plotting package.
import numpy as np #Array package.
import os #This package is necessary to handle saving our data.

#################################################ESTABLISHING VARIABLES
#PDE related
mx = 500*m #meters - model length
my = 100*m #meters - model width
ndx = 500 # mesh steps in x direction 
ndy = 1 # mesh steps in y direction
boundloc = mx/2 # location of boundary between two blocks
q=0.*Celsius #our heat source temperature is now zero
Tref=2273.*Celsius # Kelvin -the starting temperature of our RHS Block
rho = 2750*kg/m**3 #kg/m^{3} density of granite
cp = 790.*J/(kg*K) #j/Kg.K thermal capacity
rhocp = rho*cp	#DENSITY * SPECIFIC HEAT
eta=0.  # RADIATION CONDITION
kappa=2.2*W/m/K #watts/m.K thermal conductivity

#Script/Iteration Related
t=0. #our start time, usually zero
tend=10*yr #the time we want to end the simulation in years
outputs = 400 # number of time steps required.
h=(tend-t)/outputs #size of time step
#user warning statement
print "Expected Number of Output Files is: ", outputs
print "Step size is: ", h/(24.*60*60), "days"
i=0 #loop counter 
#the folder to put our outputs in, leave blank "" for script path 
save_path="data/onedheatdiff002"
########## note this folder path must exist to work ###################

################################################ESTABLISHING PARAMETERS
#generate domain using rectangle
model = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
#extract finite points - the solution points
x=model.getX()
#create the PDE
mypde=LinearPDE(model) #assigns a domain to our PDE
mypde.setSymmetryOn() #set the fast solver on for symmetry
#define our PDE coeffs
mypde.setValue(A=kappa*kronecker(model),D=rhocp/h,d=eta,y=eta*Tref)
#establish location of boundary between two blocks
bound = x[0]-boundloc
#set initial temperature
T= 0*Tref*whereNegative(bound)+Tref*wherePositive(bound)

#convert solution points for plotting
plx = x.toListOfTuples() 
plx = np.array(plx) #convert to tuple to numpy array
plx = plx[:,0] #extract x locations

########################################################START ITERATION
while t<=tend:
	i+=1 #increment the counter
	t+=h #increment the current time
	mypde.setValue(Y=rhocp/h*T) #reset variable PDE coefficients
	T=mypde.getSolution() #find temperature solution
	#set up for plotting
	tempT = T.toListOfTuples(scalarastuple=False)
	pl.figure(1)
	pl.plot(plx,tempT)
	pl.axis([0,500,0,2500])
	pl.title("Temperature accross Interface")
	pl.savefig(os.path.join(save_path,"intpyplot%03d.png") %i)
	pl.clf()
	
# compile the *.png files to create an *.avi video that shows T change
# with time. This opperation uses linux mencoder.
os.system("mencoder mf://"+save_path+"/*.png -mf type=png:\
w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
onedheatdiff002tempT.avi")
