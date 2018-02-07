
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.pycad import *
from esys.pycad.gmsh import Design
from esys.finley import MakeDomain
from esys.escript import unitsSI as U

DIM=3
NE_D=10
DEPTH=1.*U.km
LX=5.*U.km
LY=2.*U.km
DIP=20*U.DEG
STRIKE=10*U.DEG

des=Design(dim=DIM, order=2, element_size = DEPTH/NE_D, keep_files=True)

if DIM ==2 :
  if (DEPTH*cos(DIP) > 0.9 * LX):
      raise ValueError("X edge to short (maybe DIP to large.)")

  p0=Point(DEPTH*cos(DIP),-DEPTH)
  p1=Point(LX,-DEPTH)
  p2=Point(LX,0.)
  p3=Point(0.,0.)

  l01=Line(p0,p1)
  l12=Line(p1,p2)
  l23=Line(p2,p3)
  l30=Line(p3,p0)
  s=PlaneSurface(CurveLoop(l01,l12,l23,l30))
  des.addItems(s, PropertySet("subduction",l30))

else: 

  if (DEPTH*cos(DIP) > 0.9 * LX):
      raise ValueError("X edge to short (maybe DIP to large.)")
  if (DEPTH*cos(DIP)+LY*sin(STRIKE) > 0.9 * LX):
      raise ValueError("X edge to short (maybe DIP to large.)")
  if (LY*sin(STRIKE) > 0.9 * LX):
      raise ValueError("X edge to short (maybe DIP to large.)")

 
  p0=Point(DEPTH*cos(DIP),0.,-DEPTH)
  p1=Point(LX,0.,-DEPTH)
  p2=Point(DEPTH*cos(DIP)+LY*sin(STRIKE),LY,-DEPTH)
  p3=Point(LX,LY,-DEPTH)
  p4=Point(0.,0.,0.)
  p5=Point(LX,0.,0.)
  p6=Point(LY*sin(STRIKE),LY,0.)
  p7=Point(LX,LY,0.)
  
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
  des.addItems(v, PropertySet("subduction", left))

des.setScriptFileName("sub.geo")
des.setMeshFileName("sub.msh")

dom=MakeDomain(des, useMacroElements=True)
dom.write("sub.fly")
