from __future__ import division
from __future__ import print_function
##############################################################################
#
# Copyright (c) 2009-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2009-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

############################################################FILE HEADER
# example09.py
# Antony Hallam
# Seismic Wave Equation Simulation using acceleration solution.
# 3D model with multiple layers.

#######################################################EXTERNAL MODULES
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.escript import *
from esys.finley import Rectangle
from esys.weipa import saveVTK
import os
# smoothing operator 
from esys.escript.pdetools import Projector, Locator
from esys.escript.unitsSI import *
import numpy as np

import pylab as pl
import matplotlib.cm as cm
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.finley import ReadMesh

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
    import sys
    print("This example will not run in an MPI world.")
    sys.exit(0)

#################################################ESTABLISHING VARIABLES
# where to save output data
savepath = "data/example09c"
meshpath = "data/example09n"
mkDir(savepath)
#Geometric and material property related variables.
domain=ReadMesh(os.path.join(savepath,'example09n.fly')) # create the domain
x=Solution(domain).getX()
#parameters layers 1,2,3,4 and fault
prho=np.array([2200.,2500.,3200.,4500.,5500.]) #density
pvel=np.array([1500.,2200.,3000.,3200.,5000.]) #velocity
pmu=pvel**2.*prho/4.                              #bulk modulus
plam=pvel**2.*prho/2.                             #lames constant
nlayers=4
width=300.0
rho=Scalar(0,Function(domain))
vel=Scalar(0,Function(domain))
mu=Scalar(0,Function(domain))
lam=Scalar(0,Function(domain))

print(0.5*np.sqrt(prho/(plam+2*pmu))*0.5)

for i in range(0,nlayers):
    rho.setTaggedValue('lblock%d'%i,prho[i])
    rho.setTaggedValue('rblock%d'%i,prho[i])
    vel.setTaggedValue('lblock%d'%i,pvel[i])
    vel.setTaggedValue('rblock%d'%i,pvel[i])
    mu.setTaggedValue('lblock%d'%i,pmu[i])
    mu.setTaggedValue('rblock%d'%i,pmu[i])
    lam.setTaggedValue('lblock%d'%i,plam[i])
    lam.setTaggedValue('rblock%d'%i,plam[i])
i=nlayers
rho.setTaggedValue('fault',prho[i])
vel.setTaggedValue('fault',pvel[i])
mu.setTaggedValue('fault',pmu[i])
lam.setTaggedValue('fault',plam[i])


# Time related variables.
testing=False
if testing:
    print('The testing end time is currently selected. This severely limits the number of time iterations.')
    print("Try changing testing to False for more iterations.")
    tend=0.1
else:
    tend=0.1    # end time

h=0.00001    # time step
# data recording times
rtime=0.0 # first time to record
rtime_inc=tend/750.0 # time increment to record
#Check to make sure number of time steps is not too large.
print("Time step size= ",h, "Expected number of outputs= ",tend/h)

U0=0.1 # amplitude of point source
ls=500   # length of the source
source=np.zeros(ls,'float') # source array
decay1=np.zeros(ls,'float') # decay curve one
decay2=np.zeros(ls,'float') # decay curve two
time=np.zeros(ls,'float')   # time values
g=np.log(0.01)/ls

dfeq=50 #Dominant Frequency
a = 2.0 * (np.pi * dfeq)**2.0
t0 = 5.0 / (2.0 * np.pi * dfeq)
srclength = 5. * t0
ls = int(srclength/h)
print('source length',ls)
source=np.zeros(ls,'float') # source array
ampmax=0
for it in range(0,ls):
    t = it*h
    tt = t-t0
    dum1 = np.exp(-a * tt * tt)
    source[it] = -2. * a * tt * dum1
    if (abs(source[it]) > ampmax):
        ampmax = abs(source[it])
    time[t]=t*h

# will introduce a spherical source at middle left of bottom face
xc=[150,0]

##########################################################ESTABLISH PDE
mypde=LinearPDE(domain) # create pde
mypde.setSymmetryOn() # turn symmetry on
# turn lumping on for more efficient solving
mypde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
kmat = kronecker(domain) # create the kronecker delta function of the domain
mypde.setValue(D=rho*kmat) #set the general form value D

############################################FIRST TIME STEPS AND SOURCE
# define small radius around point xc
src_length = 10; print("src_length = ",src_length)
# set initial values for first two time steps with source terms
xb=FunctionOnBoundary(domain).getX()
yx=(cos(length(xb-xc)*3.1415/src_length)+1)*whereNegative(length(xb-xc)-src_length)
stop=Scalar(0.0,FunctionOnBoundary(domain))
stop.setTaggedValue("top",1.0)
src_dir=numpy.array([0.,-1.]) # defines direction of point source as down

mypde.setValue(y=source[0]*yx*src_dir*stop) #set the source as a function on the boundary

# initial value of displacement at point source is constant (U0=0.01)
# for first two time steps
u=[0.0,0.0]*x
u_m1=u

####################################################ITERATION VARIABLES
n=0 # iteration counter
t=0 # time counter
##############################################################ITERATION
while t<tend:
    # get current stress
    g=grad(u); stress=lam*trace(g)*kmat+mu*(g+transpose(g))#*abc
    mypde.setValue(X=-stress) # set PDE values
    accel = mypde.getSolution() #get PDE solution for accelleration
    u_p1=(2.*u-u_m1)+h*h*accel #calculate displacement
    u_p1=u_p1#*abc          # apply boundary conditions
    u_m1=u; u=u_p1 # shift values by 1
    # save current displacement, acceleration and pressure
    if (t >= rtime):
        saveVTK(os.path.join(savepath,"ex09c.%05d.vtu"%n),displacement=length(u),\
                                    acceleration=length(accel),tensor=stress)
        rtime=rtime+rtime_inc #increment data save time
    # increment loop values
    t=t+h; n=n+1
    if (n < ls):
        mypde.setValue(y=source[n]*yx*src_dir*stop) #set the source as a function on the boundary
    print("time step %d, t=%s"%(n,t))
