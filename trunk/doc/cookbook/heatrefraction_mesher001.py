
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

from esys.pycad import *
from esys.pycad.gmsh import Design
from esys.finley import MakeDomain
import numpy as np
#import numpy as n
from math import *

# routine to find consecutive coordinates of a loop in pycad
def loopcoords(loop):
	# return all construction points of input
	temp = loop.getConstructionPoints()
	#create a numpy array for xyz components or construction points
	coords = np.zeros([len(temp),3],float)
	#place construction points in array
	for i in range(0,len(temp)):
		coords[i,:]=temp[i].getCoordinates()
	#return a numpy array
	return coords

# Overall Domain
p0=Point(0.0,        0.0, 0.0)
p1=Point(0.0,    -6000.0, 0.0)
p2=Point(5000.0, -6000.0, 0.0)
p3=Point(5000.0,     0.0, 0.0)

l01=Line(p0, p1)
l12=Line(p1, p2)
l23=Line(p2, p3)
l30=Line(p3, p0)

c=CurveLoop(l01, l12, l23, l30)

# Material Boundary

x=[ Point(i*100.0,-2500+1500.*cos(pi*i*100.0/2500.0+pi)) for i in range(0,51) ]
mysp = Spline(*tuple(x))
x1=Spline.getStartPoint(mysp)
x2=Spline.getEndPoint(mysp)
	
# TOP BLOCK
tbl1=Line(p0,x1)
tbl2=mysp
tbl3=Line(x2,p3)
tblockloop = CurveLoop(tbl1,tbl2,tbl3,l30)
tblock = PlaneSurface(tblockloop)

tpg = loopcoords(tblockloop)
np.savetxt("toppg",tpg,delimiter=" ")

# BOTTOM BLOCK
bbl1=Line(x1,p1)
bbl3=Line(p2,x2)
bbl4=-mysp
bblockloop = CurveLoop(bbl1,l12,bbl3,bbl4)
bblock = PlaneSurface(bblockloop)

#clockwise check
bblockloop2=CurveLoop(mysp,Line(x2,p2),Line(p2,p1),Line(p1,x1))
bpg = loopcoords(bblockloop2)
np.savetxt("botpg",bpg,delimiter=" ")

# Create a Design which can make the mesh
d=Design(dim=2, element_size=200)
# Add the trapezoid with cutout
d.addItems(PropertySet("top",tblock),PropertySet("bottom",bblock),PropertySet("linebottom",l12))
# Create the geometry, mesh and Escript domain
d.setScriptFileName("heatrefraction_mesh001.geo")

d.setMeshFileName("heatrefraction_mesh001.msh")
domain=MakeDomain(d, integrationOrder=-1, reducedIntegrationOrder=-1, optimizeLabeling=True)
# Create a file that can be read back in to python with mesh=ReadMesh(fileName)
domain.write("heatrefraction_mesh001.fly")


