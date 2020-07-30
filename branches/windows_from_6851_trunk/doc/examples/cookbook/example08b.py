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
# example08b.py
# Antony Hallam
# Seismic Wave Equation Simulation using acceleration solution.
# Extend the solution in example 08a to use absorbing boundary
# conditions.

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
    savepath = "data/example08b"
    mkDir(savepath)
    #Geometric and material property related variables.
    mx = 1000. # model lenght
    my = 1000. # model width
    ndx = 300 # steps in x direction 
    ndy = 300 # steps in y direction
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

    h=0.0001    # time step
    # data recording times
    rtime=0.0 # first time to record
    rtime_inc=tend/50.0 # time increment to record
    #Check to make sure number of time steps is not too large.
    print("Time step size= ",h, "Expected number of outputs= ",tend/h)

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
    #   source[it] = exp(-a * tt * tt)    !gaussian
        if (abs(source[it]) > ampmax):
            ampmax = abs(source[it])
        #source[t]=np.exp(g*t)*U0*np.sin(2.*np.pi*t/(0.75*ls))*(np.exp(-.1*g*t)-1)
        #decay1[t]=np.exp(g*t)
        #decay2[t]=(np.exp(-.1*g*t)-1)
        time[it]=t*h
    #tdecay=decay1*decay2*U0
    #decay1=decay1*U0; decay2=decay2*U0
    pl.clf(); 
    pl.plot(source)
    #pl.plot(time,decay1);pl.plot(time,decay2); 
    #pl.plot(time,tdecay)
    pl.savefig(os.path.join(savepath,'source.png'))

    # will introduce a spherical source at middle left of bottom face
    xc=[mx/2,0]

    ####################################################DOMAIN CONSTRUCTION
    domain=Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy,order=2) # create the domain
    x=domain.getX() # get the locations of the nodes in the domani

    ##########################################################ESTABLISH PDE
    mypde=LinearPDE(domain) # create pde
    mypde.setSymmetryOn() # turn symmetry on
    # turn lumping on for more efficient solving
    mypde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
    kmat = kronecker(domain) # create the kronecker delta function of the domain
    mypde.setValue(D=kmat*rho) #set the general form value D

    ##########################################################ESTABLISH ABC
    # Define where the boundary decay will be applied.
    bn=50.
    bleft=xstep*bn; bright=mx-(xstep*bn); bbot=my-(ystep*bn)
    # btop=ystep*bn # don't apply to force boundary!!!

    # locate these points in the domain
    left=x[0]-bleft; right=x[0]-bright; bottom=x[1]-bbot

    tgamma=0.85   # decay value for exponential function
    def calc_gamma(G,npts):
        func=np.sqrt(abs(-1.*np.log(G)/(npts**2.)))
        return func

    gleft  = calc_gamma(tgamma,bleft)
    gright = calc_gamma(tgamma,bleft)
    gbottom= calc_gamma(tgamma,ystep*bn)

    print('gamma', gleft,gright,gbottom)

    # calculate decay functions
    def abc_bfunc(gamma,loc,x,G):
        func=exp(-1.*(gamma*abs(loc-x))**2.)
        return func

    fleft=abc_bfunc(gleft,bleft,x[0],tgamma)
    fright=abc_bfunc(gright,bright,x[0],tgamma)
    fbottom=abc_bfunc(gbottom,bbot,x[1],tgamma)
    # apply these functions only where relevant
    abcleft=fleft*whereNegative(left)
    abcright=fright*wherePositive(right)
    abcbottom=fbottom*wherePositive(bottom)
    # make sure the inside of the abc is value 1
    abcleft=abcleft+whereZero(abcleft)
    abcright=abcright+whereZero(abcright)
    abcbottom=abcbottom+whereZero(abcbottom)
    # multiply the conditions together to get a smooth result
    abc=abcleft*abcright*abcbottom

    #visualise the boundary function
    #abcT=abc.toListOfTuples()
    #abcT=np.reshape(abcT,(ndx+1,ndy+1))
    #pl.clf(); pl.imshow(abcT); pl.colorbar(); 
    #pl.savefig(os.path.join(savepath,"abc.png"))


    ############################################FIRST TIME STEPS AND SOURCE
    # define small radius around point xc
    src_length = 40; print("src_length = ",src_length)
    # set initial values for first two time steps with source terms
    y=source[0]*(cos(length(x-xc)*3.1415/src_length)+1)*whereNegative(length(x-xc)-src_length)
    src_dir=np.array([0.,1.]) # defines direction of point source as down
    y=y*src_dir
    mypde.setValue(y=y) #set the source as a function on the boundary
    # turn lumping on for more efficient solving
    mypde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
    # for first two time steps
    u=[0.0,0.0]*wherePositive(x)
    u_m1=u

    ####################################################ITERATION VARIABLES
    n=0 # iteration counter
    t=0 # time counter
    ##############################################################ITERATION
    while t<tend:
        # get current stress
        g=grad(u); stress=lam*trace(g)*kmat+mu*(g+transpose(g))
        mypde.setValue(X=-stress*abc) # set PDE values
        accel = mypde.getSolution() #get PDE solution for accelleration
        u_p1=(2.*u-u_m1)+h*h*accel #calculate displacement
        u_p1=u_p1*abc       # apply boundary conditions
        u_m1=u; u=u_p1 # shift values by 1
        # save current displacement, acceleration and pressure
        if (t >= rtime):
            saveVTK(os.path.join(savepath,"ex08b.%05d.vtu"%n),displacement=length(u),\
                                        acceleration=length(accel),tensor=stress)
            rtime=rtime+rtime_inc #increment data save time
        # increment loop values
        t=t+h; n=n+1
        if (n < ls):
            y=source[n]*(cos(length(x-xc)*3.1415/src_length)+1)*whereNegative(length(x-xc)-src_length)
            y=y*src_dir; mypde.setValue(y=y) #set the source as a function on the boundary
        print("time step %d, t=%s"%(n,t))
