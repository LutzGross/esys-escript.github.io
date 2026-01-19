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


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import subprocess

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
from esys.finley import Rectangle, Brick, Merge, JoinFaces
from esys.finley import ReadGmsh
from esys.weipa import saveVTK
import os
GEO = { 2 :
"""
Mesh.MshFileVersion = 2.2;
Mesh.ElementOrder = {o};
l_X={l_X};
element_size = {element_size};
Point(10)=\{0.,0.,0.,element_size\};
Point(11)=\{l_X,0.,0.,element_size\};
Point(12)=\{l_X,1.,0.,element_size\};
Point(13)=\{0.,1.,0.,element_size\};
Line(1)=\{10,11\};
Line(12)=\{11,12\};
Line(23)=\{12,13\};
Line(30)=\{13,10\};
CurveLoop(10) = \{1,12,23,30\};
Plane_Surface(1) = \{10\};
Physical Curve(100) = \{1, 12, 23, 30\};
Physical Surface(200) = \{1\};
""",
        3 :
"""
Mesh.MshFileVersion = 2.2;
Mesh.ElementOrder = {o};
l_X={l_X};
element_size = {element_size};
Point(1000)=\{ 0.,0.,0.,element_size\};
Point(1100)=\{l_X,0.,0.,element_size\};
Point(1010)=\{0.,1.,0.,element_size\};
Point(1110)=\{l_X,1.,0.,element_size\};
Point(1001)=\{0.,0.,1.,element_size\};
Point(1101)=\{l_X,0.,1.,element_size\};
Point(1011)=\{0.,1.,1.,element_size\};
Point(1111)=\{l_X,1.,1.,element_size\};

Line(10)=\{1000,1100\};
Line(20)=\{1100,1110\};
Line(30)=\{1110,1010\};
Line(40)=\{1010,1000\};
Line(11)=\{1000,1001\};
Line(21)=\{1100,1101\};
Line(31)=\{1110,1111\};
Line(41)=\{1010,1011\};
Line(12)=\{1001,1101\};
Line(22)=\{1101,1111\};
Line(32)=\{1111,1011\};
Line(42)=\{1011,1001\};
CurveLoop(10)=\{10,20,30,40\};
CurveLoop(11)=\{12,22,32,42\};
CurveLoop(12)=\{10,21,-12,-11\}; 
CurveLoop(13)=\{30,41,-32,-31\};
CurveLoop(14)=\{11,-42,-41,40\};
CurveLoop(15)=\{-21,20,31,-22\};

PlaneSurface(1)=\{-10\};
PlaneSurface(2)=\{11\};
PlaneSurface(3)=\{12}\;
PlaneSurface(4)=\{13}\;
PlaneSurface(5)=\{14\};
PlaneSurface(6)=\{15\};
SurfaceLoop(10) = \{1,2,3,4,5,6\}; 
Volume(1) = \{10\};
Physical Surface(200) = \{1,2,3,4,5,6\}; 
Physical Volume(300) = \{1\};
"""
        }
def getMesh(NE_X, NE_Y, t,d,o,fullOrder,r,l_X):
   if t == "Hex":
       if d == 2:
         dom=Rectangle(n0=NE_X, n1=NE_Y, l0=l_X,order=o, useFullElementOrder=fullOrder,useElementsOnFace=r,optimize=True)
       else:
         Brick()
         dom=Brick(n0=NE_X, n1=NE_Y, n2=NE_Y,l0=l_X,order=o, useFullElementOrder=fullOrder,useElementsOnFace=r,optimize=True)
   else:
      geospec = GEO[d].format( {"o": o, "l_X" : l_X,  "element_size" : min(l_X/max(3,NE_X),1./max(3,NE_Y)) } )
      FF=open("tester.geo", 'w')
      FF.write(geospec)
      FF.close()
      import subprocess
      rp = subprocess.run(["gmsh", f"-{o}", "-optimize_netgen", "-o", "tester.msh", FF.name], stdout=subprocess.PIPE)
      dom = ReadGmsh("tester.msh", d, optimize=True)
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
