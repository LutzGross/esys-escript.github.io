
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
# onedheatdiff001.py
# Model temperature diffusion in an Iron bar. This is a one dimensional
# problem with a single heat source at the LHS

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
import cblib

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

#################################################ESTABLISHING VARIABLES
#Domain related.
mx = 1*m #meters - model length
my = .1*m #meters - model width
ndx = 100 # mesh steps in x direction 
ndy = 1 # mesh steps in y direction - one dimension means one element

#PDE related
q=200. * Celsius #Kelvin - our heat source temperature
Tref = 0. * Celsius #Kelvin - starting temp of iron bar
rho = 7874. *kg/m**3 #kg/m^{3} density of iron
cp = 449.*J/(kg*K) #j/Kg.K thermal capacity
rhocp = rho*cp
kappa = 80.*W/m/K #watts/m.Kthermal conductivity

#Script/Iteration Related
t=0 #our start time, usually zero
tend=5.*minute #seconds - time to end simulation
outputs = 200 # number of time steps required.
h=(tend-t)/outputs #size of time step
#user warning statement
print "Expected Number of time outputs is: ", (tend-t)/h
i=0 #loop counter
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","onedheatdiff001")

#ensure the dir exists
cblib.needdirs([save_path, os.path.join(save_path,"tempT"),\
                         os.path.join(save_path, "totT")])

################################################ESTABLISHING PARAMETERS
#generate domain using rectangle
rod = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
#extract finite points - the solution points
x=rod.getX()
#create the PDE
mypde=LinearPDE(rod) #assigns a domain to our PDE
mypde.setSymmetryOn() #set the fast solver on for symmetry
mypde.setValue(A=kappa*kronecker(rod),D=rhocp/h) #define our PDE coeffs
qH=q*whereZero(x[0]) #set heat source
T=Tref # set initial temperature

#convert solution points for plotting
plx = x.toListOfTuples() 
plx = np.array(plx) #convert to tuple to numpy array
plx = plx[:,0] #extract x locations

########################################################START ITERATION
while t<=tend:
	i+=1 #increment the counter
	t+=h #increment the current time
	mypde.setValue(Y=qH+rhocp/h*T) #set variable PDE coefficients
	T=mypde.getSolution() #get the PDE solution
	totT = rhocp*T #get the total heat solution in the system	
	
	#establish figure 1 for temperature vs x plots
	tempT = T.toListOfTuples(scalarastuple=False)
	pl.figure(1) #current figure
	pl.plot(plx,tempT) #plot solution
		#define axis extents and title
	pl.axis([0,1.0,273.14990+0.00008,0.004+273.1499])
	pl.title("Temperature accross Rod")
		#save figure to file
	if getMPIRankWorld() == 0:
		pl.savefig(os.path.join(save_path,"tempT",\
		                                  "rodpyplot%03d.png"%i))
	pl.clf() #clear figure
	
	#establish figure 2 for total temperature vs x plots and repeat
	tottempT = totT.toListOfTuples(scalarastuple=False)
	pl.figure(2)
	pl.plot(plx,tottempT)
	pl.axis([0,1.0,9.657E08,12000+9.657E08])
	pl.title("Total temperature accross Rod")
	if getMPIRankWorld() == 0:
		pl.savefig(os.path.join(save_path,"totT",\
		                                  "ttrodpyplot%03d.png"%i))
	pl.clf()

# compile the *.png files to create two *.avi videos that show T change
# with time. This opperation uses linux mencoder. For other operating 
# systems it is possible to use your favourite video compiler to
# convert image files to videos. To enable this step uncomment the
# following lines.

#os.system("mencoder mf://"+save_path+"/tempT"+"/*.png -mf type=png:\
#w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
#onedheatdiff001tempT.avi")

#os.system("mencoder mf://"+save_path+"/totT"+"/*.png -mf type=png:\
#w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
#onedheatdiff001totT.avi")
