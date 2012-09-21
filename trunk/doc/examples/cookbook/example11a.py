
##############################################################################
#
# Copyright (c) 2009-2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2009-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example10a.py
# Model of electrical potential in a half-space.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.weipa import saveVTK
from esys.finley import Rectangle # This imports the rectangle domain function from finley
import os, sys #This package is necessary to handle saving our data.
from math import pi, sqrt, sin, cos

import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.escript.pdetools import Projector
from cblib import toRegGrid


import pylab as pl #Plotting package
import numpy as np

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
    import sys
    print("This example will not run in an MPI world.")
    sys.exit(0)

#################################################ESTABLISHING VARIABLES
#Domain related.
mx = 2000*m #meters - model length
my = -1000*m #meters - model depth
ndx = 200 # mesh steps in x direction 
ndy = 100 # mesh steps in y direction - one dimension means one element
#PDE related
res=1000.0
con=1/res
cur=10.

################################################ESTABLISHING PARAMETERS
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","example11")
#ensure the dir exists
mkDir(save_path)

####################################################DOMAIN CONSTRUCTION
domain = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
x=Solution(domain).getX()

kro=kronecker(domain)
source1=[mx/4.,0]; source2=[3.*mx/4.,0]

sourceg=length(exp(-length(x-source1)/(100.)))+length(exp(-length(x-source2)/(100.)))
sourceg=sourceg/integrate(sourceg)

q=whereZero(x[1]-my)+whereZero(x[0])+whereZero(x[0]-mx)
###############################################ESCRIPT PDE CONSTRUCTION

mypde=LinearPDE(domain)
mypde.setValue(A=kro*con,Y=sourceg,q=q,r=0)
mypde.setSymmetryOn()
sol=mypde.getSolution()

# Save the output to file.
saveVTK(os.path.join(save_path,"ex11a.vtu"),source=sourceg,res_pot=sol)
