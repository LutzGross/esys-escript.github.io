from __future__ import division, print_function
##############################################################################
#
# Copyright (c) 2009-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2009-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

############################################################FILE HEADER
# example09.py
# Antony Hallam
# Seismic Wave Equation Simulation using acceleration solution.
# 3D model with multiple layers. Layercake example.

#######################################################EXTERNAL MODULES
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.escript import *
from esys.weipa import saveVTK
import os
# smoothing operator 
from esys.escript.pdetools import Projector, Locator
from esys.escript.unitsSI import *
import numpy as np

import pylab as pl
import matplotlib.cm as cm
from esys.escript.linearPDEs import LinearPDE, SolverOptions
try:
    # This imports the rectangle domain function 
    from esys.finley import Rectangle, ReadMesh
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
        import sys
        print("This example will not run in an MPI world.")
        sys.exit(0)

if HAVE_FINLEY:
    #################################################ESTABLISHING VARIABLES
    # where to save output data
    savepath = "data/example09b"
    meshpath = "data/example09m"
    mkDir(savepath)
    #Geometric and material property related variables.
    step=4.0 # the element size

    vel=1800.    #starting velocity
    rhoc=2000.   #starting density
    nlayers=9  #number of layers in layercake model.

    ####################################################TESTING SWITCH
    testing=True
    if testing:
        print('The testing end time is currently selected. This severely limits the number of time iterations.')
        print("Try changing testing to False for more iterations.")
        tend=0.001
        #Model Parameters    
        mx=40.
        my=40.
        mz=20.
        outputs=5
    else:
        tend=0.1    # end time
        #Model Parameters
        mx=100.0   #x width of model
        my=100.0   #y width of model
        mz=50.0   #depth of model
        outputs=200

    ####################################################TIME RELATED VARIABLES 
    h=0.00001    # time step
    # data recording times
    rtime=0.0 # first time to record
    rtime_inc=tend/outputs # time increment to record
    #Check to make sure number of time steps is not too large.
    print("Time step size= ",h, "Expected number of outputs= ",tend/h)

    ####################################################CREATING THE SOURCE FUNCTION
    U0=0.1 # amplitude of point source
    dfeq=50 #Dominant Frequency
    a = 2.0 * (np.pi * dfeq)**2.0
    t0 = 5.0 / (2.0 * np.pi * dfeq)
    srclength = 5. * t0

    ls = int(srclength/h)
    print('source length',ls)

    source=np.zeros(ls,'float') # source array
    decay1=np.zeros(ls,'float') # decay curve one
    decay2=np.zeros(ls,'float') # decay curve two
    time=np.zeros(ls,'float')   # time values
    g=np.log(0.01)/ls

    ampmax=0
    for it in range(0,ls):
        t = it*h
        tt = t-t0
        dum1 = np.exp(-a * tt * tt)
        source[it] = -2. * a * tt * dum1
        if (abs(source[it]) > ampmax):
            ampmax = abs(source[it])
        time[it]=t*h

    # will introduce a spherical source at middle left of bottom face
    xc=[mx/2,my/2,0]

    ####################################################DOMAIN CONSTRUCTION
    domain=ReadMesh(os.path.join(meshpath,'example09lc.fly')) # create the domain
    x=domain.getX() # get the locations of the nodes in the domain

    lam=Scalar(0,Function(domain))
    mu=Scalar(0,Function(domain))
    rho=Scalar(0,Function(domain))

    #Setting parameters for each layer in the model.
    for i in range(0,nlayers):
        rho.setTaggedValue("volume_%d"%i,rhoc+i*100.)
        lamc=(vel+i*100.)**2.*(rhoc+i*100.)/2.
        muc=(vel+i*100.)**2.*(rhoc+i*100.)/4.
        lam.setTaggedValue("volume_%d"%i,lamc)
        mu.setTaggedValue("volume_%d"%i,muc)

    ##########################################################ESTABLISH PDE
    mypde=LinearPDE(domain) # create pde
    mypde.setSymmetryOn() # turn symmetry on
    # turn lumping on for more efficient solving
    #mypde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
    kmat = kronecker(domain) # create the kronecker delta function of the domain
    mypde.setValue(D=rho*kmat) #set the general form value D

    ############################################FIRST TIME STEPS AND SOURCE
    # define small radius around point xc
    src_rad = 20; print("src radius= ",src_rad)
    # set initial values for first two time steps with source terms
    xb=FunctionOnBoundary(domain).getX()
    yx=(cos(length(xb-xc)*3.1415/src_rad)+1)*whereNegative(length(xb-xc)-src_rad)
    stop=Scalar(0.0,FunctionOnBoundary(domain))
    stop.setTaggedValue("intface_0",1.0)
    src_dir=numpy.array([0.,0.,1.0]) # defines direction of point source as down

    mypde.setValue(y=source[0]*yx*src_dir*stop) #set the source as a function on the boundary
    # initial value of displacement at point source is constant (U0=0.01)
    # for first two time steps
    u=[0.0,0.0,0.0]*wherePositive(x)
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
            saveVTK(os.path.join(savepath,"ex09b.%05d.vtu"%n),displacement=length(u),\
                                        acceleration=length(accel),tensor=stress)
            rtime=rtime+rtime_inc #increment data save time
        # increment loop values
        t=t+h; n=n+1
        if (n < ls):
            mypde.setValue(y=source[n]*yx*src_dir*stop) #set the source as a function on the boundary
        print("time step %d, t=%s"%(n,t))
