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

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example08c.py
# Create either a 2D syncline or anticline model using pycad meshing 
# tools. Wave equation solution.

#######################################################EXTERNAL MODULES
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.pycad import * #domain constructor
from esys.pycad.gmsh import Design #Finite Element meshing package
import os #file path tool
from math import * # math package
from esys.escript import *
from esys.escript.unitsSI import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.pdetools import Projector
from esys.weipa import saveVTK
from cblib import toRegGrid, subsample
import matplotlib
matplotlib.use('agg') #It's just here for automated testing

import pylab as pl #Plotting package
import numpy as np
try:
    from esys.finley import MakeDomain
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
    #set modal to 1 for a syncline or -1 for an anticline structural 
    #configuration
    modal=-1

    # the folder to put our outputs in, leave blank "" for script path - 
    # note this folder path must exist to work
    save_path= os.path.join("data","example08c") 
    mkDir(save_path)

    ################################################ESTABLISHING PARAMETERS
    #Model Parameters
    width=1000.0*m   #width of model
    depth=1000.0*m  #depth of model
    dx=5
    xstep=dx # calculate the size of delta x
    ystep=dx # calculate the size of delta y

    sspl=51 #number of discrete points in spline
    dsp=width/(sspl-1) #dx of spline steps for width
    dep_sp=500.0*m #avg depth of spline
    h_sp=250.*m #heigh of spline
    orit=-1.0 #orientation of spline 1.0=>up -1.0=>down

    vel2=1800.;   vel1=3000.
    rho2=2300.;   rho1=3100. #density
    mu2=(vel2**2.)*rho2/8.;  mu1=(vel1**2.)*rho1/8.  #bulk modulus
    lam2=mu2*6.; lam1=mu1*6. #lames constant


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
    # will introduce a spherical source at middle left of bottom face
    xc=[width/2,0]
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
        if (abs(source[it]) > ampmax):
            ampmax = abs(source[it])
        time[it]=t*h

    ####################################################DOMAIN CONSTRUCTION
    # Domain Corners
    p0=Point(0.0,      0.0, 0.0)
    p1=Point(0.0,    depth, 0.0)
    p2=Point(width, depth, 0.0)
    p3=Point(width,   0.0, 0.0)

    # Generate Material Boundary
    x=[ Point(i*dsp\
        ,dep_sp+modal*orit*h_sp*cos(pi*i*dsp/dep_sp+pi))\
         for i in range(0,sspl)\
      ]
    mysp = Spline(*tuple(x))
    # Start and end of material boundary.
    x1=mysp.getStartPoint()
    x2=mysp.getEndPoint()

    #  Create TOP BLOCK
    # lines
    tbl1=Line(p0,x1)
    tbl2=mysp
    tbl3=Line(x2,p3)
    l30=Line(p3, p0)
    # curve
    tblockloop = CurveLoop(tbl1,tbl2,tbl3,l30)
    # surface
    tblock = PlaneSurface(tblockloop)
    # Create BOTTOM BLOCK
    # lines
    bbl1=Line(x1,p1)
    bbl3=Line(p2,x2)
    bbl4=-mysp
    l12=Line(p1, p2)
    # curve
    bblockloop = CurveLoop(bbl1,l12,bbl3,bbl4)

    # surface
    bblock = PlaneSurface(bblockloop)

    #clockwise check as splines must be set as polygons in the point order
    #they were created. Otherwise get a line across plot.
    bblockloop2=CurveLoop(mysp,Line(x2,p2),Line(p2,p1),Line(p1,x1))

    ################################################CREATE MESH FOR ESCRIPT
    # Create a Design which can make the mesh
    d=Design(dim=2, element_size=dx, order=2)
    # Add the subdomains and flux boundaries.
    d.addItems(PropertySet("top",tblock),PropertySet("bottom",bblock),\
                                         PropertySet("linetop",l30))
    # Create the geometry, mesh and Escript domain
    d.setScriptFileName(os.path.join(save_path,"example08c.geo"))
    d.setMeshFileName(os.path.join(save_path,"example08c.msh"))
    domain=MakeDomain(d, optimizeLabeling=True)
    x=domain.getX()
    print("Domain has been generated ...")

    lam=Scalar(0,Function(domain))
    lam.setTaggedValue("top",lam1)
    lam.setTaggedValue("bottom",lam2)
    mu=Scalar(0,Function(domain))
    mu.setTaggedValue("top",mu1)
    mu.setTaggedValue("bottom",mu2)
    rho=Scalar(0,Function(domain))
    rho.setTaggedValue("top",rho1)
    rho.setTaggedValue("bottom",rho2)

    ##########################################################ESTABLISH PDE
    mypde=LinearPDE(domain) # create pde
    mypde.setSymmetryOn() # turn symmetry on
    # turn lumping on for more efficient solving
    #mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
    kmat = kronecker(domain) # create the kronecker delta function of the domain
    mypde.setValue(D=rho*kmat) #set the general form value D

    ##########################################################ESTABLISH ABC
    # Define where the boundary decay will be applied.
    bn=20.
    bleft=xstep*bn; bright=width-(xstep*bn); bbot=depth-(ystep*bn)
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

    ############################################FIRST TIME STEPS AND SOURCE
    # define small radius around point xc
    src_length = 40; print("src_length = ",src_length)
    # set initial values for first two time steps with source terms
    xb=FunctionOnBoundary(domain).getX()
    y=source[0]*(cos(length(x-xc)*3.1415/src_length)+1)*whereNegative(length(xb-src_length))
    src_dir=numpy.array([0.,1.]) # defines direction of point source as down
    y=y*src_dir
    mypde.setValue(y=y) #set the source as a function on the boundary
    # initial value of displacement at point source is constant (U0=0.01)
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
            saveVTK(os.path.join(save_path,"ex08c.%05d.vtu"%n),\
                        vector_displacement=u,displacement=length(u),\
                        vector_acceleration=accel,acceleration=length(accel),\
                        tensor=stress)
            rtime=rtime+rtime_inc #increment data save time
        # increment loop values
        t=t+h; n=n+1
        if (n < ls):
            y=source[n]*(cos(length(x-xc)*3.1415/src_length)+1)*whereNegative(length(x-xc)-src_length)
            y=y*src_dir; mypde.setValue(y=y) #set the source as a function on the boundary
        print("time step %d, t=%s"%(n,t))
