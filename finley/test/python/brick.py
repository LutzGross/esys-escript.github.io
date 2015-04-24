
##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2003-2015 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.pycad import *
from esys.pycad.gmsh import Design
from esys.finley import MakeDomain


p0=Point(0.,0.,0.)
p1=Point(1.,0.,0.)
p2=Point(0.,1.,0.)
p3=Point(1.,1.,0.)
p4=Point(0.,0.,1.)
p5=Point(1.,0.,1.)
p6=Point(0.,1.,1.)
p7=Point(1.,1.,1.)

l01=Line(p0,p1)
l13=Line(p1,p3)
l32=Line(p3,p2)
l20=Line(p2,p0)

l45=Line(p4,p5)
l57=Line(p5,p7)
l76=Line(p7,p6)
l64=Line(p6,p4)

l15=Line(p1,p5)
l40=Line(p4,p0)
l37=Line(p3,p7)
l62=Line(p6,p2)

bottom=PlaneSurface(-CurveLoop(l01,l13,l32,l20))
top=PlaneSurface(CurveLoop(l45,l57,l76,l64))
front=PlaneSurface(CurveLoop(l01,l15,-l45,l40))
back=PlaneSurface(CurveLoop(l32,-l62,-l76,-l37))
left=PlaneSurface(CurveLoop(-l40,-l64,l62,l20))
right=PlaneSurface(CurveLoop(-l15,l13,l37,-l57))
v=Volume(SurfaceLoop(top,bottom,front,back,left,right))
des=Design(dim=3, order=2, element_size = 1, keep_files=True)
des.addItems(PropertySet("V",v), PropertySet("bottom",bottom))
des.setScriptFileName("dom.geo")
des.setMeshFileName("dom.msh")

dom=MakeDomain(des)
dom.write("brick.fly")



