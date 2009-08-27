
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

from esys.pycad import * #domain constructor
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.finley import MakeDomain #Converter for escript
import os #file path tool
import numpy as np #numerial python for arrays
from math import * # math package
#used to construct polygons for plotting from pycad
from cblib import getLoopCoords 

#set modal to 1 for a syncline or -1 for an anticline structural 
#configuration	
modal=1

# the folder to put our outputs in, leave blank "" for script path - 
# note this folder path must exist to work
save_path = "data/heatrefrac001" 

#Model Parameters
width=5000.0   #width of model
depth=-6000.0  #depth of model
sspl=51 #number of discrete points in spline
dsp=width/(sspl-1) #dx of spline steps for width
dep_sp=2500.0 #avg depth of spline
h_sp=1500.0 #heigh of spline
orit=-1.0 #orientation of spline 1.0=>up -1.0=>down

# Overall Domain
p0=Point(0.0,      0.0, 0.0)
p1=Point(0.0,    depth, 0.0)
p2=Point(width, depth, 0.0)
p3=Point(width,   0.0, 0.0)

l01=Line(p0, p1)
l12=Line(p1, p2)
l23=Line(p2, p3)
l30=Line(p3, p0)

c=CurveLoop(l01, l12, l23, l30)

# Material Boundary
x=[ Point(i*dsp\
    ,-dep_sp+modal*orit*h_sp*cos(pi*i*dsp/dep_sp+pi))\
     for i in range(0,sspl)\
  ]
     
mysp = Spline(*tuple(x))
x1=Spline.getStartPoint(mysp)
x2=Spline.getEndPoint(mysp)
	
# TOP BLOCK
tbl1=Line(p0,x1)
tbl2=mysp
tbl3=Line(x2,p3)
tblockloop = CurveLoop(tbl1,tbl2,tbl3,l30)
tblock = PlaneSurface(tblockloop)

tpg = getLoopCoords(tblockloop)
np.savetxt(os.path.join(save_path,"toppg"),tpg,delimiter=" ")

# BOTTOM BLOCK
bbl1=Line(x1,p1)
bbl3=Line(p2,x2)
bbl4=-mysp
bblockloop = CurveLoop(bbl1,l12,bbl3,bbl4)
bblock = PlaneSurface(bblockloop)

#clockwise check as splines must be set as polygons in the point order
#they were created. Otherwise get a line across plot.
bblockloop2=CurveLoop(mysp,Line(x2,p2),Line(p2,p1),Line(p1,x1))
bpg = getLoopCoords(bblockloop2)
np.savetxt(os.path.join(save_path,"botpg"),bpg,delimiter=" ")

# Create a Design which can make the mesh
d=Design(dim=2, element_size=200)
# Add the subdomains and flux boundaries.
d.addItems(PropertySet("top",tblock),PropertySet("bottom",bblock),PropertySet("linebottom",l12))
# Create the geometry, mesh and Escript domain
d.setScriptFileName(os.path.join(save_path,"heatrefraction_mesh001.geo"))

d.setMeshFileName(os.path.join(save_path,"heatrefraction_mesh001.msh"))
domain=MakeDomain(d, integrationOrder=-1, reducedIntegrationOrder=-1, optimizeLabeling=True)
# Create a file that can be read back in to python with mesh=ReadMesh(fileName)
domain.write(os.path.join(save_path,"heatrefraction_mesh001.fly"))


