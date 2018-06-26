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
# Model of gravitational Potential.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.weipa import saveVTK
import os, sys #This package is necessary to handle saving our data.

try:
    # This imports the rectangle domain function 
    from esys.finley import ReadMesh
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
    print("This example will not run in an MPI world.")
    sys.exit(0)

if HAVE_FINLEY:
    #################################################ESTABLISHING VARIABLES
    G=6.67300*10E-11

    ################################################ESTABLISHING PARAMETERS
    #the folder to put our outputs in, leave blank "" for script path 
    save_path= os.path.join("data","example10")
    #ensure the dir exists
    mkDir(save_path)

    ####################################################DOMAIN CONSTRUCTION
    domain=ReadMesh(os.path.join(save_path,'fault.fly')) # create the domain
    x=Solution(domain).getX()
    rho=Scalar(0,Function(domain))
    rho.setTaggedValue("xx",500.)
    rho.setTaggedValue("limestone",0.0)
    rho.setTaggedValue("fault",1200.)

    kro=kronecker(domain)

    q=whereZero(x[2])#-sup(x[2]))
    ###############################################ESCRIPT PDE CONSTRUCTION

    mypde=LinearPDE(domain)
    mypde.setValue(A=kro,Y=4.*3.1415*G*rho,q=q,r=0)
    sol=mypde.getSolution()
    saveVTK(os.path.join(save_path,"ex10c.vtu"),\
            grav_pot=sol,\
            g_field=-grad(sol),\
            g_fieldz=-grad(sol)*[0,0,1],\
            gz=length(-grad(sol)*[0,0,1]),\
            fault=rho)
