
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# You can shorten the execution time by reducing variable tend from 60 to 0.5

# Importing all the necessary modules required.
from esys.escript import *
from esys.finley import Rectangle
import sys
import os
from cblib import *

# Establish a save path.
savepath = "data/wavesolver2d001nwtest"
# Creating a directory automatically to store the output data.
if not os.path.isdir("data"):
   os.mkdir("data")
if not os.path.isdir(savepath):
   os.mkdir(savepath)


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
print "Time step size= ",h, "Expected number of outputs= ",tend/h
proceeder = raw_input("Is this ok?(y/n)")
#Exit if user thinks too many outputs.
if proceeder == "n":
   sys.exit()

U0=0.01 # amplitude of point source
#  spherical source at middle of bottom face

xc=[300,200]

mydomain=Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
wavesolver2d(mydomain,h,tend,lam,mu,rho,U0,xc,savepath)

