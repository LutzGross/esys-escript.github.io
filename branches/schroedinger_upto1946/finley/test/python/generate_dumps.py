########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"
"""
this script generates dump files of various data type (scalar, vector, tensor)
on various function spaces. 

The meshes are read from the MESH_DIRECTORY. filenames need to start with "mesh_" and 
have the extension ".fly". Besides the mesh dump files (extension ".nc") and the
coresponding data files with extensions (_<fs_name>_s.nc, _<fs_name>_v.nc and _<fs_name>_t.nc where
<fs_name> is the name of the function space)

This script can run under MPI.

"""

MESH_DIRECTORY="./tmp_meshes"
from esys.escript import *
from esys.finley import ReadMesh
import os

for root, dirs, files in os.walk(MESH_DIRECTORY, topdown=False):
   for name in files: 
       f=name.split(".")
       if f[0].startswith("mesh_") and f[-1]=="fly":
          print "start reading ",os.path.join(MESH_DIRECTORY,name)
          dom=ReadMesh(os.path.join(MESH_DIRECTORY,name),optimize=True)
          dom.dump(os.path.join(MESH_DIRECTORY,f[0]+".nc"))
          for fs_name in ["ContinuousFunction", "Solution", "Function", "FunctionOnBoundary", "FunctionOnContactZero", "FunctionOnContactOne", 
                          "ReducedContinuousFunction", "ReducedSolution", "ReducedFunction", "ReducedFunctionOnBoundary", "ReducedFunctionOnContactZero", "ReducedFunctionOnContactOne"]:
             if fs_name == "ContinuousFunction":
                 fs= ContinuousFunction(dom)
             if fs_name == "Solution":
                 fs= Solution(dom)
             if fs_name == "Function":
                 fs= Function(dom)
             if fs_name == "FunctionOnBoundary":
                 fs= FunctionOnBoundary(dom)
             if fs_name == "FunctionOnContactZero":
                 fs == FunctionOnContactZero(dom)
             if fs_name == "FunctionOnContactOne":
                 fs = FunctionOnContactOne(dom)
             if fs_name == "ReducedContinuousFunction":
                 fs= ReducedContinuousFunction(dom)
             if fs_name == "ReducedSolution":
                 fs= ReducedSolution(dom)
             if fs_name == "ReducedFunction":
                 fs= ReducedFunction(dom)
             if fs_name == "ReducedFunctionOnBoundary":
                 fs= ReducedFunctionOnBoundary(dom)
             if fs_name == "ReducedFunctionOnContactZero":
                 fs == ReducedFunctionOnContactZero(dom)
             if fs_name == "ReducedFunctionOnContactOne":
                 fs = ReducedFunctionOnContactOne(dom)
             for type in [ "s", "v", "t" ]:
                 data_file=os.path.join(MESH_DIRECTORY,f[0]+"_"+fs_name+"_"+type+".nc")
                 x=fs.getX()
                 print "\t data file ",data_file
                 if type == "t":
                    n=normalize(x)/length(x-Vector(0.5,fs))
                    d=outer(n,x)
                 elif type == "v":
                    d=normalize(x)/length(x-Vector(0.5,fs))
                 else:
                    d=length(x-Vector(0.5,fs))
                 d.dump(data_file)
