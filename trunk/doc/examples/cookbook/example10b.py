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
# example10a.py
# Model of gravitational Potential for a Gravity Well.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.weipa import saveVTK # This imports the VTK file saver from weipa
import os, sys #This package is necessary to handle saving our data.

from esys.escript.pdetools import Projector, Locator
import numpy as np

try:
    # This imports the rectangle domain function 
    from esys.finley import Brick
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
    my = 500*m #meters - model width
    mz = -4000*m
    ndx = 15 # mesh steps in x direction 
    ndy = 15 # mesh steps in y direction - one dimension means one element
    ndz = 10
    #PDE related
    rho=100.0
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
    domain.setX(interpolate(x, ContinuousFunction(domain)))
    mask=wherePositive(100-length(x-rholoc))
    rho=rho*mask
    kro=kronecker(domain)

    mass=rho*vol(domain)
    ipot=FunctionOnBoundary(domain)
    xb=ipot.getX()

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

    ################################################MODEL SIZE SAMPLING
    sampler=[]
    for i in range(-250,250,1):
        sampler.append([i,0,250])

    sample=[] # array to hold values
    rec=Locator(domain,sampler) #location to record
    psol=rec.getValue(sol)
    np.savetxt(os.path.join(save_path,"example10b_%04d.asc"%mx),psol)

