
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
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.dudley import MakeDomain #Converter for escript
from esys.escript import mkDir, getMPISizeWorld
import os
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

# make sure path exists 
save_path= os.path.join("data","example09") 
mkDir(save_path)

################################################ESTABLISHING PARAMETERS
#Model Parameters
xwidth=2000.0   #x width of model
ywidth=2000.0   #y width of model
depth=500.0   #depth of model
intf=depth/2.   #Depth of the interface.

####################################################DOMAIN CONSTRUCTION
# Domain Corners
p0=Point(0.0,    0.0,      0.0)
p1=Point(xwidth, 0.0,      0.0)
p2=Point(xwidth, ywidth,   0.0)
p3=Point(0.0,    ywidth,   0.0)
p4=Point(0.0,    ywidth, depth)
p5=Point(0.0,    0.0,    depth)
p6=Point(xwidth, 0.0,    depth)
p7=Point(xwidth, ywidth, depth)
# Join corners in anti-clockwise manner.
l01=Line(p0, p1)
l12=Line(p1, p2)
l23=Line(p2, p3)
l30=Line(p3, p0)
l56=Line(p5, p6)
l67=Line(p6, p7)
l74=Line(p7, p4)
l45=Line(p4, p5)

# Join line segments to create domain boundaries and then surfaces
ctop=CurveLoop(l01, l12, l23, l30);     stop=PlaneSurface(ctop)
cbot=CurveLoop(-l67, -l56, -l45, -l74); sbot=PlaneSurface(cbot)

# for each side
ip0=Point(0.0,    0.0,      intf)
ip1=Point(xwidth, 0.0,      intf)
ip2=Point(xwidth, ywidth,   intf)
ip3=Point(0.0,    ywidth,   intf)

linte_ar=[]; #lines for vertical edges
linhe_ar=[]; #lines for horizontal edges
linte_ar.append(Line(p0,ip0))
linte_ar.append(Line(ip0,p5))
linte_ar.append(Line(p1,ip1))
linte_ar.append(Line(ip1,p6))
linte_ar.append(Line(p2,ip2))
linte_ar.append(Line(ip2,p7))
linte_ar.append(Line(p3,ip3))
linte_ar.append(Line(ip3,p4))

linhe_ar.append(Line(ip0,ip1))
linhe_ar.append(Line(ip1,ip2))
linhe_ar.append(Line(ip2,ip3))
linhe_ar.append(Line(ip3,ip0))

cintfa_ar=[]; cintfb_ar=[] #curveloops for above and below interface on sides
cintfa_ar.append(CurveLoop(linte_ar[0],linhe_ar[0],-linte_ar[2],-l01))
cintfa_ar.append(CurveLoop(linte_ar[2],linhe_ar[1],-linte_ar[4],-l12))
cintfa_ar.append(CurveLoop(linte_ar[4],linhe_ar[2],-linte_ar[6],-l23))
cintfa_ar.append(CurveLoop(linte_ar[6],linhe_ar[3],-linte_ar[0],-l30))

cintfb_ar.append(CurveLoop(linte_ar[1],l56,-linte_ar[3],-linhe_ar[0]))
cintfb_ar.append(CurveLoop(linte_ar[3],l67,-linte_ar[5],-linhe_ar[1]))
cintfb_ar.append(CurveLoop(linte_ar[5],l74,-linte_ar[7],-linhe_ar[2]))
cintfb_ar.append(CurveLoop(linte_ar[7],l45,-linte_ar[1],-linhe_ar[3]))

sintfa_ar=[PlaneSurface(cintfa_ar[i]) for i in range(0,4)]
sintfb_ar=[PlaneSurface(cintfb_ar[i]) for i in range(0,4)]

sintf=PlaneSurface(CurveLoop(*tuple(linhe_ar)))

vintfa=Volume(SurfaceLoop(stop,-sintf,*tuple(sintfa_ar)))
vintfb=Volume(SurfaceLoop(sbot,sintf,*tuple(sintfb_ar)))

# Create the volume.
#sloop=SurfaceLoop(stop,sbot,*tuple(sintfa_ar+sintfb_ar))
#model=Volume(sloop)


#############################################EXPORTING MESH FOR ESCRIPT
# Create a Design which can make the mesh
d=Design(dim=3, element_size=200.0)

d.addItems(PropertySet('vintfa',vintfa))
d.addItems(PropertySet('vintfb',vintfb))
d.addItems(sintf)

d.setScriptFileName(os.path.join(save_path,"example09n.geo"))
d.setMeshFileName(os.path.join(save_path,"example09n.msh"))
#
#  
#
domain=MakeDomain(d)
# Create a file that can be read back in to python with
# mesh=ReadMesh(fileName)
domain.write(os.path.join(save_path,"example09n.fly"))

from esys.pycad import layer_cake
intfaces=[10,30,50,55,80,100,200,250,400,500]

cmplx_domain=layer_cake.LayerCake(xwidth,ywidth,intfaces,200.0)
cmplx_domain.setScriptFileName(os.path.join(save_path,"example09lc.geo"))
cmplx_domain.setMeshFileName(os.path.join(save_path,"example09lc.msh"))
dcmplx=MakeDomain(cmplx_domain)
dcmplx.write(os.path.join(save_path,"example09lc.fly"))




