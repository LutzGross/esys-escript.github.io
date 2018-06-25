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
# example03b.py
# Model temperature diffusion between a granite intrusion and sandstone 
# country rock. This is a two dimensional problem with the granite as a
# heat source. It creates vtk files.
# 
#  This program is MPI safe.

#######################################################EXTERNAL MODULES
#To solve the problem it is necessary to import the modules we require.
#This imports everything from the escript library
from esys.escript import *
# This defines the LinearPDE module as LinearPDE
from esys.escript.linearPDEs import LinearPDE 
# This imports the VTK file saver function
from esys.weipa import saveVTK
# A useful unit handling package which will make sure all our units
# match up in the equations under SI.
from esys.escript.unitsSI import *
import os
try:
    # This imports the rectangle domain function 
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False

if HAVE_FINLEY:
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

    ########################################################START ITERATION
    while t<=tend:
          i+=1 #counter
          t+=h #current time
          mypde.setValue(Y=qH+T*rhocp/h)
          T=mypde.getSolution()
          saveVTK(os.path.join(save_path,"data.%03d.vtu"%i), T=T)
          print("time step %s at t=%e days completed."%(i,t/day))

    # use 
    #
    #  cd data/example03
    #  mayavi2 -d data.001.vtu -m Surface
    #
    # to visualize the results (mayavi2 must be installed on your system).
    #
