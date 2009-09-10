
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

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
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
# heatrefraction_mesher002.py
# Create a 3 block fault and overburden style model in 2d using pycad
# meshing tools.

#######################################################EXTERNAL MODULES
from esys.pycad import *
from esys.pycad.gmsh import Design
from esys.finley import MakeDomain
import os
import numpy as np
from math import *
from cblib import needdirs, getLoopCoords

#################################################ESTABLISHING VARIABLES
# where to put output files
save_path= os.path.join("data","heatrefrac002")
needdirs([save_path])

################################################ESTABLISHING PARAMETERS
# Overall Domain
p0=Point(0.0,        0.0, 0.0)
p1=Point(0.0,    -6000.0, 0.0)
p2=Point(5000.0, -6000.0, 0.0)
p3=Point(5000.0,     0.0, 0.0)

####################################################DOMAIN CONSTRUCTION
l01=Line(p0, p1)
l12=Line(p1, p2)
l23=Line(p2, p3)
l30=Line(p3, p0)

c=CurveLoop(l01, l12, l23, l30)

# Generate Material Boundary
p4=Point(0.0,    -2400.0, 0.0)
p5=Point(2000.0, -2400.0, 0.0)
p6=Point(3000.0, -6000.0, 0.0)
p7=Point(5000.0, -2400.0, 0.0)
	
# Create TOP BLOCK
tbl1=Line(p0,p4)
tbl2=Line(p4,p5)
tbl3=Line(p5,p7)
tbl4=Line(p7,p3)
tblockloop = CurveLoop(tbl1,tbl2,tbl3,tbl4,l30)
tblock = PlaneSurface(tblockloop)
# python plotting polygon
tpg = loopcoords(tblockloop)
np.savetxt(os.path.join(save_path,"toppg"),tpg,delimiter=" ")

# Create BOTTOM BLOCK LEFT
bbll1=Line(p4,p1)
bbll2=Line(p1,p6)
bbll3=Line(p6,p5)
bbll4=-tbl2
bblockloopl = CurveLoop(bbll1,bbll2,bbll3,bbll4)
bblockl = PlaneSurface(bblockloopl)

#clockwise check
#bblockloopl2=CurveLoop(mysp,Line(x2,p2),Line(p2,p1),Line(p1,x1))
bpg = loopcoords(bblockloopl)
np.savetxt(os.path.join(save_path,"botpgl"),bpg,delimiter=" ")

# Create BOTTOM BLOCK RIGHT
bbrl1=Line(p6,p2)
bbrl2=Line(p2,p7)
bbrl3=-tbl3
bbrl4=-bbll3
bblockloopr = CurveLoop(bbrl1,bbrl2,bbrl3,bbrl4)
bblockr = PlaneSurface(bblockloopr)

#clockwise check
#bblockloopr2=CurveLoop(mysp,Line(x2,p2),Line(p2,p1),Line(p1,x1))
bpg = loopcoords(bblockloopr)
np.savetxt(os.path.join(save_path,"botpgr"),bpg,delimiter=" ")

#############################################EXPORTING MESH FOR ESCRIPT
# Create a Design which can make the mesh
d=Design(dim=2, element_size=200)
# Add the subdomains and flux boundaries.
d.addItems(PropertySet("top",tblock),\
                       PropertySet("bottomleft",bblockl),\
                       PropertySet("bottomright",bblockr),\
                       PropertySet("linebottom",bbll21, bbrl1))

# Create the geometry, mesh and Escript domain
d.setScriptFileName(os.path.join(save_path,\
                      "heatrefraction_mesh003.geo"))

d.setMeshFileName(os.path.join(save_path,"heatrefraction_mesh003.msh"))
domain=MakeDomain(d, integrationOrder=-1,\
                     reducedIntegrationOrder=-1,\
                     optimizeLabeling=True)
# Create a file that can be read back in to python with 
# mesh=ReadMesh(fileName)
domain.write(os.path.join(save_path,"heatrefraction_mesh003.fly"))

