# -*- coding: utf-8 -*-

########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Layer Cake

This script uses the pycad module to build a rectangular layered
system of prisms for modelling.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Antony Hallam antony.hallam@uqconnect.edu.au"

from esys.pycad import * #domain constructor
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.finley import MakeDomain #Converter for escript
from esys.escript.unitsSI import *
import os
from math import *
import pylab as pl
import numpy as np



def buildFreeSurface(xwidth,ywidth):
    # Build a free surface for start of layer cake model.
    # Layer Corners
    corner_points=[]
    corner_points.append(Point(0.0,    0.0,      0.0))
    corner_points.append(Point(xwidth, 0.0,      0.0))
    corner_points.append(Point(xwidth, ywidth,   0.0))
    corner_points.append(Point(0.0,    ywidth,   0.0))
    corner_points.append(corner_points[0])
    hor_lines=[]
    for i in range(0,4): #top of layer loop
        hor_lines.append(Line(corner_points[i],corner_points[i+1]))
    #free surface
    free_surf = PlaneSurface(CurveLoop(*tuple(hor_lines[0:4])))
    return free_surf,hor_lines,corner_points

def buildLayer(xwidth,ywidth,depth,lay_surf,hor_lines,corner_points):
    # Layer Corners
    corner_points.append(Point(0.0,    0.0,    depth))
    corner_points.append(Point(xwidth, 0.0,    depth))
    corner_points.append(Point(xwidth, ywidth, depth))
    corner_points.append(Point(0.0,    ywidth, depth))
    corner_points.append(corner_points[5])
    print corner_points

    # Join corners in anti-clockwise manner.
#    for i in range(0,4): #top of layer loop
#        hor_lines.append(Line(corner_points[i],corner_points[i+1])
    for i in range(0,4): #bottom of layer loop
        hor_lines.append(Line(corner_points[5+i],corner_points[6+i]))
    # Join corners vertically.
    ver_lines=[]
    for i in range(0,5):
        ver_lines.append(Line(corner_points[i],corner_points[i+5]))
        print 'corners',i,i+5

    #bottom of layer
    lay_surf=[PlaneSurface(CurveLoop(-hor_lines[3],-hor_lines[2],\
                                     -hor_lines[1],-hor_lines[1]))]
    lay_surf.append(PlaneSurface(CurveLoop(*tuple(hor_lines[4:8]))))
#    lay_subot=PlaneSurface(-CurveLoop(*tuple(hor_lines[4:8])))
    #sides of layer
    for i in range(0,4):
        lay_surf.append(PlaneSurface(CurveLoop(\
                        -ver_lines[i],\
                        -hor_lines[i+4],\
                        ver_lines[i+1],\
                        hor_lines[i]          )))
    lay_vol=Volume(-SurfaceLoop(*tuple(lay_surf)))
    #lay_vol=lay_surf
    return lay_vol,lay_surf[1],hor_lines[4:8],corner_points[5:10]

def LayerCake(xwidth,ywidth,depths,ele_size,fname,save_path=""):
    #get number of layers
    ndepths=len(depths)
    fsuf,fsurl,fsurp=buildFreeSurface(xwidth,ywidth)
    domain=Design(dim=3,element_size=ele_size)
    #domain.addItems(PropertySet('free_surf',fsuf))
    #build each layer sequentially
    tsuf=fsuf; tsurl=fsurl; tsurp=fsurp
    for i in range(0,ndepths):
        tvol,tsuf,tsurl,tsurp=buildLayer(xwidth,ywidth,depths[i],\
                                     tsuf,tsurl,tsurp)
        #domain.addItems(*tuple(tvol))
        #domain.addItems(PropertySet('intface_%d'%(i+1),tsuf))
        domain.addItems(PropertySet('volume_%d'%i,tvol))
    # Add the subdomains and flux boundaries.

    domain.setScriptFileName(os.path.join(save_path,fname+".geo"))
    domain.setMeshFileName(os.path.join(save_path,fname+".msh"))
    findomain=MakeDomain(domain) #  make the finley domain:
    # Create a file that can be read back in to python with
    #mesh=ReadMesh(fileName)
    #findomain.write(os.path.join(save_path,fname+".fly"))    

LayerCake(100.0,100.0,[100.],10.,'testlc')

