
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
# Model of gravitational Potential for a Gravity Well.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.finley import Brick # This imports the rectangle domain function from finley
from esys.weipa import saveVTK # This imports the VTK file saver from weipa
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
mz = -500*m
ndx = 50 # mesh steps in x direction 
ndy = 50 # mesh steps in y direction - one dimension means one element
ndz = 50
#PDE related
rho=1000.0
rholoc=[0,0,-0]
G=6.67300*10E-11

################################################ESTABLISHING PARAMETERS
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","example10")
#ensure the dir exists
mkDir(save_path)

####################################################DOMAIN CONSTRUCTION
domain = Brick(l0=mx,l1=my,n0=ndx, n1=ndy,l2=mz,n2=ndz)
x=Solution(domain).getX()
x=x-[mx/2,my/2,mz/2]
domain.setX(x)
mask=wherePositive(100-length(x-rholoc))
rho=rho*mask
kro=kronecker(domain)

q=whereZero(x[2]-inf(x[2]))
###############################################ESCRIPT PDE CONSTRUCTION

mypde=LinearPDE(domain)
mypde.setValue(A=kro,Y=4.*3.1415*G*rho,q=q,r=0)
mypde.setSymmetryOn()
sol=mypde.getSolution()
saveVTK(os.path.join(save_path,"ex10b.vtu"),\
        grav_pot=sol,\
        g_field=-grad(sol),\
        g_fieldz=-grad(sol)*[0,0,1],\
        gz=length(-grad(sol)*[0,0,1]))

