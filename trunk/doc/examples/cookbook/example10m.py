
##############################################################################
#
# Copyright (c) 2009-2013 by University of Queensland
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

__copyright__="""Copyright (c) 2009-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""


############################################################FILE HEADER
# example10m.py
# Create a simple 2D mesh, which is optimised for cells close to the 
# source. Larger elements are used to decrease computational requirements
# and to properly fullfil the boundary conditions.
#
#######################################################EXTERNAL MODULES
from esys.pycad import * #domain constructor
from esys.pycad.extras import layer_cake
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.finley import MakeDomain #Converter for escript
from esys.escript import mkDir, getMPISizeWorld
import os
import subprocess as sp
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
        import sys
        print("This example will not run in an MPI world.")
        sys.exit(0)

# make sure path exists 
save_path= os.path.join("data","example10m") 
mkDir(save_path)

################################################BIG DOMAIN
#ESTABLISHING PARAMETERS
width=10000.    #width of model
depth=10000.    #depth of model
bele_size=500.  #big element size
#DOMAIN CONSTRUCTION
p0=Point(0.0,    0.0)
p1=Point(width, 0.0)
p2=Point(width, depth)
p3=Point(0.0,   depth)
# Join corners in anti-clockwise manner.
l01=Line(p0, p1)
l12=Line(p1, p2)
l23=Line(p2, p3)
l30=Line(p3, p0)

cbig=CurveLoop(l01,l12,l23,l30)

################################################SMALL DOMAIN
#ESTABLISHING PARAMETERS
xwidth=1000.0   #x width of model
zdepth=1000.0   #y width of model
sele_size=10.   #small element size
#TRANSFORM
xshift=width/2-xwidth/2
zshift=depth/2-zdepth/2
#DOMAIN CONSTRUCTION
p4=Point(xshift, zshift)
p5=Point(xwidth+xshift, zshift)
p6=Point(xwidth+xshift, zdepth+zshift)
p7=Point(xshift,    zdepth+zshift)
# Join corners in anti-clockwise manner.
l45=Line(p4, p5)
l56=Line(p5, p6)
l67=Line(p6, p7)
l74=Line(p7, p4)

csmall=CurveLoop(l45,l56,l67,l74)

ssmall=PlaneSurface(csmall)
sbig=PlaneSurface(cbig,holes=[csmall])

#############################################EXPORTING MESH FOR ESCRIPT
# Design the geometry for the big mesh.
d1=Design(dim=2, element_size=bele_size, order=1)
d1.addItems(sbig)
d1.addItems(PropertySet(l01,l12,l23,l30))
d1.setScriptFileName(os.path.join(save_path,"example10m_big.geo"))
MakeDomain(d1)

# Design the geometry for the small mesh.
d2=Design(dim=2, element_size=sele_size, order=1)
d2.addItems(ssmall)
d2.setScriptFileName(os.path.join(save_path,"example10m_small.geo"))
MakeDomain(d2)

# Join the two meshes using Gmsh and then apply a 2D meshing algorithm.
# The small mesh must come before the big mesh in the merging call!!@!!@!
sp.call("gmsh -2 "+
        os.path.join(save_path,"example10m_small.geo")+" "+
        os.path.join(save_path,"example10m_big.geo")+" -o "+
        os.path.join(save_path,"example10m.msh"),shell=True)
