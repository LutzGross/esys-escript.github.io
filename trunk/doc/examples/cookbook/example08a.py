from __future__ import division, print_function
##############################################################################
#
# Copyright (c) 2009-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2009-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

############################################################FILE HEADER
# example08a.py
# Antony Hallam
# Seismic Wave Equation Simulation using acceleration solution.

#######################################################EXTERNAL MODULES
from esys.escript import *
from esys.weipa import saveVTK
import sys
import os
# smoothing operator 
from esys.escript.pdetools import Projector, Locator
from esys.escript.unitsSI import *
import numpy as np
from esys.escript.linearPDEs import LinearPDE, SolverOptions
try:
    # This imports the rectangle domain function 
    from esys.finley import Rectangle
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
    savepath = "data/example08a"
    mkDir(savepath)
    #Geometric and material property related variables.
    mx = 1000. # model lenght
    my = -1000. # model width
    ndx = 500 # steps in x direction 
    ndy = 500 # steps in y direction
    xstep=mx/ndx # calculate the size of delta x
    ystep=abs(my/ndy) # calculate the size of delta y
    lam=3.462e9 #lames constant
    mu=3.462e9  #bulk modulus
    rho=1154.   #density
    # Time related variables.
    testing=True
    if testing:
            print('The testing end time is currently selected. This severely limits the number of time iterations.')
            print("Try changing testing to False for more iterations.")
            tend=0.001
    else:
            tend=0.5    # end time

    h=0.0005     # time step
    # data recording times
    rtime=0.0 # first time to record
    rtime_inc=tend/20.0 # time increment to record
    #Check to make sure number of time steps is not too large.
    print("Time step size= ",h, "Expected number of outputs= ",tend/h)

    U0=0.01 # amplitude of point source
    # will introduce a spherical source at middle left of bottom face
    xc=[mx/2,0]

    ####################################################DOMAIN CONSTRUCTION
    domain=Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy) # create the domain
    x=domain.getX() # get the locations of the nodes in the domani

    ##########################################################ESTABLISH PDE
    mypde=LinearPDE(domain) # create pde
    mypde.setSymmetryOn() # turn symmetry on
    # turn lumping on for more efficient solving
    mypde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
    kmat = kronecker(domain) # create the kronecker delta function of the domain
    mypde.setValue(D=kmat*rho) #set the general form value D

    ############################################FIRST TIME STEPS AND SOURCE
    # define small radius around point xc
    src_length = 20; print("src_length = ",src_length)
    # set initial values for first two time steps with source terms
    y=U0*(cos(length(x-xc)*3.1415/src_length)+1)*whereNegative(length(x-xc)-src_length)
    src_dir=np.array([0.,-1.]) # defines direction of point source as down
    y=y*src_dir
    mypde.setValue(y=y) #set the source as a function on the boundary
    # initial value of displacement at point source is constant (U0=0.01)
    # for first two time steps
    u=[0.0,0.0]*whereNegative(x)
    u_m1=u

    ####################################################ITERATION VARIABLES
    n=0 # iteration counter
    t=0 # time counter
    ##############################################################ITERATION
    while t<tend:
        # get current stress
        g=grad(u); stress=lam*trace(g)*kmat+mu*(g+transpose(g))
        mypde.setValue(X=-stress) # set PDE values
        accel = mypde.getSolution() #get PDE solution for accelleration
        u_p1=(2.*u-u_m1)+h*h*accel #calculate displacement
        u_m1=u; u=u_p1 # shift values by 1
        # save current displacement, acceleration and pressure
        if (t >= rtime):
            saveVTK(os.path.join(savepath,"ex08a.%05d.vtu"%n),displacement=length(u),\
                                        acceleration=length(accel),tensor=stress)
            rtime=rtime+rtime_inc #increment data save time
        # increment loop values
        t=t+h; n=n+1
        print("time step %d, t=%s"%(n,t))
