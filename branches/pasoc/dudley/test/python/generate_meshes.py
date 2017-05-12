##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"
"""
this script generates a variety of meshes in the fly format on a unit square with a 
fixed number of elements NE (approximatively).

This script will not run under MPI.

generated file names are:

   mesh_D<d>_o<O>_T<t>_Contacts<c>_Rich<r>.fly

wehere <d>= spatial dimension 2,3
       <o>= order 1,2, 2F
       <t>= element type = hex, tet
       <c>= Yes if contact elements are present, No otherwise.
       <r>= YEs is rich surface elements are used, No otherwise

any item with x[d]<CUT will be tagged with 1 the rest will be marked 100.
"""

MESH_DIRECTORY="./tmp_meshes"
NE=10
CUT=0.5

from esys.escript import *
from esys.dudley import Rectangle, Brick, Merge, JoinFaces
from esys.pycad import Point, Line,PlaneSurface, CurveLoop, Volume,SurfaceLoop
from esys.pycad.gmsh import Design
from esys.dudley import MakeDomain
import os

def getMesh(NE_X, NE_Y, t,d,o,fullOrder,r,l_X):
   if t == "Hex":
       if d == 2:
         dom=Rectangle(n0=NE_X, n1=NE_Y, l0=l_X,order=o, useFullElementOrder=fullOrder,useElementsOnFace=r,optimize=True)
       else:
         Brick()
         dom=Brick(n0=NE_X, n1=NE_Y, n2=NE_Y,l0=l_X,order=o, useFullElementOrder=fullOrder,useElementsOnFace=r,optimize=True)
   else:
       des=Design(dim=d, order=o, element_size =min(l_X/max(3,NE_X),1./max(3,NE_Y)), keep_files=True)
       des.setScriptFileName("tester.geo")
       if d == 2:
          p0=Point(0.,0.)
          p1=Point(l_X,0.)
          p2=Point(l_X,1.)
          p3=Point(0.,1.)
          l01=Line(p0,p1)
          l12=Line(p1,p2)
          l23=Line(p2,p3)
          l30=Line(p3,p0)
          s=PlaneSurface(CurveLoop(l01,l12,l23,l30))
          des.addItems(s,l01,l12,l23,l30)
       else: 
          p000=Point( 0.,0.,0.)
          p100=Point(l_X,0.,0.)
          p010=Point(0.,1.,0.)
          p110=Point(l_X,1.,0.)
          p001=Point(0.,0.,1.)
          p101=Point(l_X,0.,1.)
          p011=Point(0.,1.,1.)
          p111=Point(l_X,1.,1.)

          l10=Line(p000,p100)
          l20=Line(p100,p110)
          l30=Line(p110,p010)
          l40=Line(p010,p000)

          l11=Line(p000,p001)
          l21=Line(p100,p101)
          l31=Line(p110,p111)
          l41=Line(p010,p011)

          l12=Line(p001,p101)
          l22=Line(p101,p111)
          l32=Line(p111,p011)
          l42=Line(p011,p001)

          bottom=PlaneSurface(-CurveLoop(l10,l20,l30,l40))
          top=PlaneSurface(CurveLoop(l12,l22,l32,l42))

          front=PlaneSurface(CurveLoop(l10,l21,-l12,-l11)) 
          back=PlaneSurface(CurveLoop(l30,l41,-l32,-l31))

          left=PlaneSurface(CurveLoop(l11,-l42,-l41,l40))
          right=PlaneSurface(CurveLoop(-l21,l20,l31,-l22))

          vol=Volume(SurfaceLoop(bottom,top,front,back,left,right))
          des.addItems(vol)

       dom=MakeDomain(des)
   return dom
    
# test if out put dir exists:
if not os.path.isdir(MESH_DIRECTORY): os.mkdir(MESH_DIRECTORY)

for d in [2,3]:
   for  o in ["1", "2", "2F"]:
      for t in [ "Hex", "Tet" ]:
         for c in ["Yes", "No"]:
            for r in ["Yes", "No"]:
               filename="mesh_D%s_o%s_T%s_Contacts%s_Rich%s.fly"%(d,o,t,c,r)
               # certain cases are excluded:
               if  ( o == "2F" and t == "Tet" ) or \
                   ( t == "Tet" and r == "Yes" ) or \
                   ( o == "2F" and r == "Yes" )  :
                 pass
               else:
                  print("generate file ",filename)
                  if c == "Yes":
                     NE_X=int(NE**(1./d)/2+0.5)
                     NE_Y=int(NE**(1./d)+0.5)
                  else:
                     NE_X=int(NE**(1./d)+0.5)
                     NE_Y=NE_X
                  print(filename)
                  print("generating ",NE_X, NE_Y)
                  if o == "2":
                     o2=2
                     full=False
                  elif o == "2F":
                     o2=2
                     full=True
                  else:
                     o2=1
                     full=False
                  if c == "Yes":
                     dom1=getMesh(NE_X, NE_Y,t,d,o2,full,r=="Yes",0.5)
                     dom2=getMesh(NE_X, NE_Y,t,d,o2,full,r=="Yes",0.5)
                     x=dom2.getX().copy()
                     x[0]=1.-x[0]
                     dom2.setX(x)
                     dom=JoinFaces([dom1,dom2])
                  else:
                     dom=getMesh(NE_X, NE_Y,t,d,o2,full,r=="Yes",1.)
                  # set tags:
                  for fs in [ContinuousFunction(dom), Function(dom), FunctionOnBoundary(dom), FunctionOnContactOne(dom)]:
                       m=whereNegative(fs.getX()[d-1]-CUT)
                       fs.getDomain().setTagMap('tag1',1)
                       fs.setTags('tag1',m)
                       fs.setTags(100,1-m)
                  dom.write(os.path.join(MESH_DIRECTORY,filename))
                  # saveVTK(os.path.join(MESH_DIRECTORY,filename+".vtu"),dom)
