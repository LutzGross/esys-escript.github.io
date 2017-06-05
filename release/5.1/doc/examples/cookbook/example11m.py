from __future__ import division, print_function
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
# example09m.py
# Create a simple 3D model for use in example09. This is the low res
# mesh for illustration purposes only.
#
#######################################################EXTERNAL MODULES
from esys.pycad import * #domain constructor
from esys.pycad.extras import layer_cake
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.escript import mkDir, getMPISizeWorld
import os
try:
    from esys.finley import MakeDomain
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
    # make sure path exists 
    save_path= os.path.join("data","example11") 
    mkDir(save_path)

    ################################################ESTABLISHING PARAMETERS
    #Model Parameters
    xwidth=500.0   #x width of model
    ywidth=500.0   #y width of model
    depth=250.0   #depth of model
    element_size=5.0

    intfaces=[50,100,200,250]

    #Specify the domain.
    domaindes=Design(dim=3,element_size=element_size,order=1)
    cmplx_domain=layer_cake(domaindes,xwidth,ywidth,intfaces)
    cmplx_domain.setScriptFileName(os.path.join(save_path,"example11lc.geo"))
    cmplx_domain.setMeshFileName(os.path.join(save_path,"example11lc.msh"))
    dcmplx=MakeDomain(cmplx_domain)
    dcmplx.write(os.path.join(save_path,"example11lc.fly"))




