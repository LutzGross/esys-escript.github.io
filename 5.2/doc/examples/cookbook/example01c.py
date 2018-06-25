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
# example01c.py
# Model temperature diffusion between two granite blocks of unequal 
# initial temperature. Solve for the spatial distribution of temperature.

# To solve the problem it is necessary to import the modules we require.
#For interactive use, you can comment out the next two lines
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
import pylab as pl #Plotting package.
import numpy as np #Array package.
import os, sys #This package is necessary to handle saving our data.
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
    #Domain related.
    mx = 500*m #meters - model length
    my = 100*m #meters - model width
    ndx = 100 # mesh steps in x direction 
    ndy = 1 # mesh steps in y direction - one dimension means one element
    boundloc = mx/2 # location of boundary between the two blocks
    #PDE related
    rho = 2750. *kg/m**3 #kg/m{3} density of iron
    cp = 790.*J/(kg*K) # J/Kg.K thermal capacity
    rhocp = rho*cp
    kappa = 2.2*W/m/K # watts/m.Kthermal conductivity
    qH=0 * J/(sec*m**3) # J/(sec.m{3}) no heat source
    T1=20 * Celsius # initial temperature at Block 1
    T2=2273. * Celsius # base temperature at Block 2

    ################################################ESTABLISHING PARAMETERS
    t=0 * day  # our start time, usually zero
    tend=50 * yr # - time to end simulation
    outputs = 200 # number of time steps required.
    h=(tend-t)/outputs #size of time step
    #user warning statement
    print("Expected Number of time outputs is: ", (tend-t)/h)
    i=0 #loop counter
    #the folder to put our outputs in, leave blank "" for script path 
    save_path= os.path.join("data","example01")
    #ensure the dir exists
    mkDir(save_path, os.path.join(save_path,"tempT"))

    ####################################################DOMAIN CONSTRUCTION
    blocks = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)

    ###############################################ESCRIPT PDE CONSTRUCTION
    #... open PDE and set coefficients ...
    mypde=LinearPDE(blocks)
    mypde.setSymmetryOn()
    A=zeros((2,2))
    A[0,0]=kappa
    mypde.setValue(A=A,D=rhocp/h)
    # ... set initial temperature ....
    x=Solution(blocks).getX()
    T= T1*whereNegative(x[0]-boundloc)+T2*(1-whereNegative(x[0]-boundloc))

    # ... open a collector for the time marks and corresponding total energy
    t_list=[]
    E_list=[]
    # ... convert solution points for plotting
    plx = x.toListOfTuples() 
    plx = np.array(plx) #convert to tuple to numpy array
    plx = plx[:,0] #extract x locations
    ########################################################START ITERATION
    while t<tend:
          i+=1
          t+=h
          mypde.setValue(Y=qH+rhocp/h*T)
          T=mypde.getSolution()
          totE=integrate(rhocp*T)
          print("time step %s at t=%e days completed. total energy = %e."%(i,t/day,totE))
          t_list.append(t)
          E_list.append(totE)

          #establish figure 1 for temperature vs x plots
          tempT = T.toListOfTuples()
          pl.figure(1) #current figure
          pl.plot(plx,tempT) #plot solution
          # add title
          pl.axis([0,mx,T1*.9,T2*1.1])
          pl.title("Temperature across blocks at time %d days"%(t/day))
          pl.ylabel('Temperature (K)')
          pl.xlabel("Length (m)")
          #save figure to file
          pl.savefig(os.path.join(save_path,"tempT", "blockspyplot%03d.png"%i))
          pl.clf() #clear figure
          
    ###############################################################PLOTTING
    # plot the total energy over time:
    pl.figure(2)
    pl.plot(t_list,E_list)
    pl.title("Total Energy")
    pl.axis([0,max(t_list),0,max(E_list)*1.1])
    pl.ylabel('Energy (W)')
    pl.xlabel('Time (s)')
    pl.savefig(os.path.join(save_path,"totE_ex01c.png"))
    pl.clf()

    ###########################################################MAKE A MOVIE
    # compile the *.png files to create a*.avi video that show T change
    # with time. This opperation uses linux mencoder. For other operating 
    # systems it may be possible to use your favourite video compiler to
    # convert image files to videos. To enable this step uncomment the
    # following lines.

    #os.system("mencoder mf://"+save_path+"/tempT"+"/*.png -mf type=png:\
    #w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
    #example01tempT.avi")
