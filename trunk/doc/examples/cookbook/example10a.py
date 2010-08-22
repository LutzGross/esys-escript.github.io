
########################################################
#
# Copyright (c) 2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example10a.py
# Model of gravitational Potential.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.escript.linearPDEs import Poisson # This defines LinearPDE as LinearPDE
from esys.finley import Rectangle # This imports the rectangle domain function from finley
#For interactive use, you can comment out the next two lines
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
import pylab as pl #Plotting package.
import numpy as np #Array package.
import os, sys #This package is necessary to handle saving our data.

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
    import sys
    print "This example will not run in an MPI world."
    sys.exit(0)

#################################################ESTABLISHING VARIABLES
#Domain related.
mx = 500*m #meters - model length
my = 500*m #meters - model width
ndx = 500 # mesh steps in x direction 
ndy = 500 # mesh steps in y direction - one dimension means one element
#PDE related
rho=100.0
rholoc=[250,250]
G=6.67300*10E-11

################################################ESTABLISHING PARAMETERS
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","example10a")
#ensure the dir exists
mkDir(save_path)

####################################################DOMAIN CONSTRUCTION
domain = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
x=Solution(domain).getX()
#mask=whereZero(x[0]-rholoc[0])*whereZero(x[1]-rholoc[1])
mask=wherePositive(100-length(x-rholoc))
#rho=Scalar(rho,ContinuousFunction(domain))
rho=rho*mask
y=Scalar(0.0,FunctionOnBoundary(domain))
kro=kronecker(domain)

q=whereZero(x[0])+\
  whereZero(x[1])+\
  whereZero(x[0]-sup(x[0]))+\
  whereZero(x[1]-sup(x[1]))
###############################################ESCRIPT PDE CONSTRUCTION

mypde=LinearPDE(domain)
mypde.setValue(A=kro,Y=4.*3.1415*G*rho,q=q,r=0.)
#mypde.setSymmetryOn()
#mypde=Poisson(domain)
#mypde.setValue(f=rho*4.*3.1415*G,q=q)
sol=mypde.getSolution()
saveVTK("ex10a.vtu",field_strength=sol)#,field=-grad(sol))
