
########################################################
#
# Copyright (c) 2009-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009-2010 by University of Queensland
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
# example09m.py
# Create a simple 3D model for use in example09. This is the low res
# mesh for illustration purposes only.
#
#######################################################EXTERNAL MODULES
from esys.pycad import * #domain constructor
from esys.pycad.extras import layer_cake
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.finley import MakeDomain #Converter for escript
from esys.escript import mkDir, getMPISizeWorld
import os
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

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




