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
from __future__ import division, print_function

__copyright__="""Copyright (c) 2009-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

# You can shorten the execution time by reducing variable tend from 60 to 0.5

# Importing all the necessary modules required.
from esys.escript import *
from esys.finley import Rectangle
import sys
import os
from cblib1 import wavesolver2d

# Establish a save path.
savepath = "data/wavesolver2d001nwtest"
mkDir(savepath)


#Geometric and material property related variables.
mx = 1000 # model lenght
my = 200 # model width
ndx = 200 # steps in x direction 
ndy = 40 # steps in y direction
lam=3.462e9 #lames constant
mu=3.462e9  #bulk modulus
rho=1154.   #density
# Time related variables.
tend=0.5    #end time
#calculating )the timestep
h=(1./5.)*sqrt(rho/(lam+2*mu))*(mx/ndx)
#Check to make sure number of time steps is not too large.
print("Time step size= ",h, "Expected number of outputs= ",tend/h)

#uncomment the following lines to give the user a chance to stop
#proceeder = raw_input("Is this ok?(y/n)")
#Exit if user thinks too many outputs.
#if proceeder == "n":
#   sys.exit()

U0=0.01 # amplitude of point source
#  spherical source at middle of bottom face

xc=[300,200]

mydomain=Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
wavesolver2d(mydomain,h,tend,lam,mu,rho,U0,xc,savepath)

