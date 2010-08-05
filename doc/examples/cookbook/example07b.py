
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

############################################################FILE HEADER
# example07b.py
# Antony Hallam
# Acoustic Wave Equation Simulation using acceleration solution
# and lumping.

#######################################################EXTERNAL MODULES
from esys.escript import *
from esys.dudley import Rectangle
import sys
import os
# smoothing operator 
from esys.escript.pdetools import Projector, Locator
from esys.escript.unitsSI import *
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from esys.escript.linearPDEs import LinearPDE

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

#################################################ESTABLISHING VARIABLES
# where to save output data
savepath = "data/example07b"
mkDir(savepath) #make sure savepath exists
#Geometric and material property related variables.
mx = 1000. # model lenght
my = 1000. # model width
ndx = 200 # steps in x direction 
ndy = 200 # steps in y direction
xstep=mx/ndx # calculate the size of delta x
ystep=my/ndy # calculate the size of delta y

c=380.0*m/sec # velocity of sound in air
csq=c*c #square of c
# Time related variables.
tend=1.5    # end time
h=0.005     # time step
# data recording times
rtime=0.0 # first time to record
rtime_inc=tend/20.0 # time increment to record
#Check to make sure number of time steps is not too large.
print "Time step size= ",h, "Expected number of outputs= ",tend/h

U0=0.005 # amplitude of point source
# want a spherical source in the middle of area
xc=[500,500] # with reference to mx,my this is the source location

####################################################DOMAIN CONSTRUCTION
mydomain=Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy) #create the domain
x=mydomain.getX() #get the node locations of the domain

##########################################################ESTABLISH PDE
mypde=LinearPDE(mydomain) # create pde
# turn lumping on for more efficient solving
mypde.getSolverOptions().setSolverMethod(mypde.getSolverOptions().LUMPING)
mypde.setSymmetryOn() # turn symmetry on
mypde.setValue(D=1.) # set the value of D in the general form to 1.

############################################FIRST TIME STEPS AND SOURCE
# define small radius around point xc
src_radius = 25.
print "src_radius = ",src_radius
# set initial values for first two time steps with source terms
u=U0*(cos(length(x-xc)*3.1415/src_radius)+1)*whereNegative(length(x-xc)-src_radius)
u_m1=u
#plot source shape
cut_loc=[] #where the cross section of the source along x will be
src_cut=[] #where the cross section of the source will be
# create locations for source cross section
for i in range(ndx/2-ndx/10,ndx/2+ndx/10):
    cut_loc.append(xstep*i)
    src_cut.append([xstep*i,xc[1]])
# locate the nearest nodes to the points in src_cut
src=Locator(mydomain,src_cut)
src_cut=src.getValue(u) #retrieve the values from the nodes
# plot the x locations vs value and save the figure
pl.plot(cut_loc,src_cut)
pl.axis([xc[0]-src_radius*3,xc[0]+src_radius*3,0.,2.*U0])
pl.savefig(os.path.join(savepath,"source_line.png"))

###########################SAVING THE VALUE AT A LOC FOR EACH TIME STEP
u_rec0=[] # array to hold values
rec=Locator(mydomain,[250.,250.]) #location to record
u_rec=rec.getValue(u); u_rec0.append(u_rec) #get the first two time steps

####################################################ITERATION VARIABLES
n=0 # iteration counter
t=0 # time counter
##############################################################ITERATION
while t<tend:
    g=grad(u); pres=csq*g # get current pressure
    mypde.setValue(X=-pres) # set values in pde
    accel = mypde.getSolution() # get new acceleration
    u_p1=(2.*u-u_m1)+h*h*accel # calculate the displacement for the next time step
    u_m1=u; u=u_p1 # shift values back one time step for next iteration
    # save current displacement, acceleration and pressure
    if (t >= rtime):
        saveVTK(os.path.join(savepath,"ex07b.%i.vtu"%n),displacement=length(u),\
                                    acceleration=length(accel),tensor=pres)
        rtime=rtime+rtime_inc #increment data save time
    u_rec0.append(rec.getValue(u)) #location specific recording
    # increment loop values
    t=t+h; n=n+1
    print n,"-th time step t ",t

# save location specific recording to file
pl.savetxt(os.path.join(savepath,'u_rec.asc'),u_rec0)
