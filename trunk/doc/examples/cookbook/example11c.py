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

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example11c.py
# Model of gravitational Potential.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
import matplotlib
matplotlib.use('agg') #It's just here for automated testing

import os, sys #This package is necessary to handle saving our data.
from math import pi, sqrt, sin, cos
import pylab as pl #Plotting package
import numpy as np
from esys.escript import * # This imports everything from the escript library
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.escript.pdetools import Projector
from esys.escript.unitsSI import * 
from esys.weipa import saveVTK

from cblib import toRegGrid

try:
    from esys.finley import ReadMesh,Rectangle
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
    my = 500*m #meters - model depth
    mz = 250
    ndx = 200 # mesh steps in x direction 
    ndy = 200 # mesh steps in y direction - one dimension means one element
    #PDE related
    res1=500.0
    res2=10.0
    res3=1000.0
    res4=10000.0
    #con=1/res
    cur=1000.

    ################################################ESTABLISHING PARAMETERS
    #the folder to put our outputs in, leave blank "" for script path 
    save_path= os.path.join("data","example11")
    #ensure the dir exists
    mkDir(save_path)

    ####################################################DOMAIN CONSTRUCTION
    domain=ReadMesh(os.path.join(save_path,'example11lc.fly')) # create the domain

    res=Scalar(0,Function(domain))
    res.setTaggedValue("volume_0",res1)
    res.setTaggedValue("volume_1",res2)
    res.setTaggedValue("volume_2",res3)
    res.setTaggedValue("volume_3",res4)
    con=1/res
    x=Solution(domain).getX()

    kro=kronecker(domain)
    source1=[3.*mx/8.,my/2,0]; source2=[5.*mx/8.,my/2,0]

    c1=length(exp(-length(x-source1)/(10.))); c1=c1/integrate(c1)
    c2=-length(exp(-length(x-source2)/(10.))); c2=c2/integrate(c2)
    sourceg=cur*(c1-c2)

    q=whereZero(x[1]-my)+whereZero(x[1])+whereZero(x[0])+whereZero(x[0]-mx)+whereZero(x[2]-mz)
    ###############################################ESCRIPT PDE CONSTRUCTION

    mypde=LinearPDE(domain)
    mypde.setValue(A=con*kro,Y=sourceg,q=q,r=0)
    #mypde.setSymmetryOn()
    sol=mypde.getSolution()

    res.expand()

    # Save the output to file.
    saveVTK(os.path.join(save_path,"ex11c.vtu"),\
            source=sourceg,\
            res_pot=sol,\
            res=res,\
            curden=-con*grad(sol),\
            abscd=length(-con*grad(sol)),\
            efield=-grad(sol))

