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
from __future__ import division, print_function

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
# example03a.py
# Model temperature diffusion between a granite intrusion and sandstone 
# country rock. This is a two dimensional problem with the granite as a
# heat source.

#######################################################EXTERNAL MODULES
#To solve the problem it is necessary to import the modules we require.
#For interactive use, you can comment out the next two lines
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
#This imports everything from the escript library
from esys.escript import *
# This defines the LinearPDE module as LinearPDE
from esys.escript.linearPDEs import LinearPDE 
# A useful unit handling package which will make sure all our units
# match up in the equations under SI.
from esys.escript.unitsSI import *
import pylab as pl #Plotting package.
import numpy as np #Array package.
import os #This package is necessary to handle saving our data.
from cblib import toXYTuple, HAVE_NATGRID
try:
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

if not HAVE_NATGRID:
    print("This example requires that natgrid is available to matplotlib")

if HAVE_FINLEY and HAVE_NATGRID:
    #################################################ESTABLISHING VARIABLES
    #PDE related
    mx = 600*m #meters - model length
    my = 600*m #meters - model width
    ndx = 150 #mesh steps in x direction 
    ndy = 150 #mesh steps in y direction
    r = 200*m #meters - radius of intrusion
    ic = [300*m, 0] #centre of intrusion (meters)
    qH=0.*J/(sec*m**3) #our heat source temperature is now zero

    ## Intrusion Variables - Granite
    Ti=2273.*Celsius # Kelvin -the starting temperature of our RHS Block
    rhoi = 2750*kg/m**3 #kg/m^{3} density of granite
    cpi = 790.*J/(kg*K) #j/Kg.K thermal capacity
    rhocpi = rhoi*cpi   #DENSITY * SPECIFIC HEAT
    kappai=2.2*W/m/K #watts/m.K thermal conductivity
    ## Country Rock Variables - Sandstone
    Tc = 473*Celsius # Kelvin #the starting temperature of our country rock
    rhoc = 2000*kg/m**3 #kg/m^{3} density
    cpc = 920.*J/(kg*K) #j/kg.k specific heat
    rhocpc = rhoc*cpc #DENSITY * SPECIFIC HEAT
    kappac = 1.9*W/m/K #watts/m.K thermal conductivity

    #Script/Iteration Related
    t=0. #our start time, usually zero
    tend=200.* yr #the time we want to end the simulation
    outputs = 200 # number of time steps required.
    h=(tend-t)/outputs #size of time step
    #user warning
    print("Expected Number of Output Files is: ", outputs)
    print("Step size is: ", h/day, "days")
    i=0 #loop counter 
    #the folder to put our outputs in, leave blank "" for script path 
    save_path= os.path.join("data","example03")
    mkDir(save_path)
    ########## note this folder path must exist to work ###################

    ################################################ESTABLISHING PARAMETERS
    #generate domain using rectangle
    model = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
    #extract finite points - the solution points
    #create the PDE
    mypde=LinearPDE(model) #assigns a domain to our PDE
    mypde.setSymmetryOn() #set the fast solver on for symmetry
    #establish location of boundary between two materials
    x=Function(model).getX()
    bound = length(x-ic)-r #where the boundary will be located
    kappa = kappai*whereNegative(bound)+kappac*(1-whereNegative(bound))
    rhocp = rhocpi*whereNegative(bound)+rhocpc*(1-whereNegative(bound))
    #define our PDE coeffs
    mypde.setValue(A=kappa*kronecker(model),D=rhocp/h)
    #set initial temperature (make sure we use the right sample points)
    x=Solution(model).getX()
    bound = length(x-ic)-r #where the boundary will be located
    T= Ti*whereNegative(bound)+Tc*(1-whereNegative(bound))

    # rearrage mymesh to suit solution function space for contouring      
    coordX, coordY = toXYTuple(T.getFunctionSpace().getX())
    # create regular grid
    xi = np.linspace(0.0,mx,75)
    yi = np.linspace(0.0,my, 75)

    ########################################################START ITERATION
    while t<=tend:
          i+=1 #counter
          t+=h #current time
          mypde.setValue(Y=qH+T*rhocp/h)
          T=mypde.getSolution()
          tempT = T.toListOfTuples()
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
          pl.savefig(os.path.join(save_path,"temp%03d.png"%i))
          pl.clf()            
          print("time step %s at t=%e days completed."%(i,t/day))

    #########################################################CREATE A MOVIE
    # compile the *.png files to create an *.avi video that shows T change
    # with time. This opperation uses linux mencoder.
    os.system("mencoder mf://"+save_path+"/*.png -mf type=png:\
w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
example03tempT.avi")
