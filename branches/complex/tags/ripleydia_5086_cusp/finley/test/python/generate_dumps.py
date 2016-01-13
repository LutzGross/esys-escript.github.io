##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"
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
          print("Reading "+os.path.join(MESH_DIRECTORY,name))
          dom=ReadMesh(os.path.join(MESH_DIRECTORY,name),optimize=True)
          for fs_name in ["ContinuousFunction", "Solution", "Function", "FunctionOnBoundary", "FunctionOnContactZero", "FunctionOnContactOne", 
                          "ReducedContinuousFunction", "ReducedSolution", "ReducedFunction", "ReducedFunctionOnBoundary", "ReducedFunctionOnContactZero", "ReducedFunctionOnContactOne"]:
             if fs_name == "ContinuousFunction":
                 fs = ContinuousFunction(dom)
             if fs_name == "Solution":
                 fs = Solution(dom)
             if fs_name == "Function":
                 fs = Function(dom)
             if fs_name == "FunctionOnBoundary":
                 fs = FunctionOnBoundary(dom)
             if fs_name == "FunctionOnContactZero":
                 fs = FunctionOnContactZero(dom)
             if fs_name == "FunctionOnContactOne":
                 fs = FunctionOnContactOne(dom)
             if fs_name == "ReducedContinuousFunction":
                 fs = ReducedContinuousFunction(dom)
             if fs_name == "ReducedSolution":
                 fs = ReducedSolution(dom)
             if fs_name == "ReducedFunction":
                 fs = ReducedFunction(dom)
             if fs_name == "ReducedFunctionOnBoundary":
                 fs = ReducedFunctionOnBoundary(dom)
             if fs_name == "ReducedFunctionOnContactZero":
                 fs = ReducedFunctionOnContactZero(dom)
             if fs_name == "ReducedFunctionOnContactOne":
                 fs = ReducedFunctionOnContactOne(dom)
             x=fs.getX()
             v=normalize(x)/clip(length(x-Vector(0.5,fs)), 1.e-8)
             t=outer(v,x)
             s=length(x-Vector(0.5,fs))
             datasetName=f[0]+"_"+fs_name
             try:
                 saveESD(datasetName, MESH_DIRECTORY, s=s, v=v, t=t)
             except:
                 print("Could not save ESD file "+datasetName)

