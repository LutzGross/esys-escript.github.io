# -*- coding: utf-8 -*-

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

def buildFreeSurface(xwidth,ywidth):
    '''
    Build a free surface to start the layer cake model.
        This surface is planar.

    Parameters:
        xwidth :: width in x direction in meters.\
        ywidth :: width in y direction in meters.\
    '''

    # Layer Corners
    corner_points=[]
    corner_points.append(Point(0.0,    0.0,      0.0))
    corner_points.append(Point(xwidth, 0.0,      0.0))
    corner_points.append(Point(xwidth, ywidth,   0.0))
    corner_points.append(Point(0.0,    ywidth,   0.0))
    corner_points.append(corner_points[0]) #repeated point for line looping

    # Edges of Free Surface 
    hor_lines=[]
    for i in range(0,4): # loop through four sides
        hor_lines.append(Line(corner_points[i],corner_points[i+1]))

    # Create Free Surface
    free_surf = PlaneSurface(CurveLoop(*tuple(hor_lines[0:4])))

    # Return Surface and primative arrays.
    return free_surf,hor_lines,corner_points

def buildLayer(xwidth,ywidth,depth,lay_surf,hor_lines,corner_points):
    '''
    Builds a boxlike volume and returns primatives for layered model
    construction.

    Parameters:
        xwidth :: width in x direction in meters.\
        ywidth :: width in y direction in meters.\
        depth  :: depth to bottom of layer in meters.\
        lay_surf :: surface at top of layer (see buildFreeSurf)\
        hor_lines :: lines of lay_surf\
        corner_points :: points of hor_lines\
    '''
    # Layer Corners
    corner_points.append(Point(0.0,    0.0,    depth))
    corner_points.append(Point(xwidth, 0.0,    depth))
    corner_points.append(Point(xwidth, ywidth, depth))
    corner_points.append(Point(0.0,    ywidth, depth))
    corner_points.append(corner_points[5]) #repeated point for line looping

    # Build the bottom surface edges.
    for i in range(0,4): # loop through four edges
        hor_lines.append(Line(corner_points[5+i],corner_points[6+i]))

    # Join corners vertically.
    ver_lines=[]
    for i in range(0,4): # loop through four corners
        ver_lines.append(Line(corner_points[i],corner_points[i+5]))
    ver_lines.append(ver_lines[0]) #repeated edge for surface looping   

    # Build surface array.
    lay_surf=[-lay_surf] #Negative of top surface
    # Bottom Surface
    lay_surf.append(PlaneSurface(CurveLoop(*tuple(hor_lines[4:8]))))
    for i in range(0,4): # loop through four sides
        lay_surf.append(PlaneSurface(CurveLoop(\
                        -ver_lines[i],-hor_lines[i+4],\
                        ver_lines[i+1],hor_lines[i] )))

    # Build Layer Volume
    lay_vol=Volume(-SurfaceLoop(*tuple(lay_surf)))
    
    # Return new volume, and primatives for next volume layer.    
    return lay_vol,-lay_surf[1],hor_lines[4:8],corner_points[5:10]


def layer_cake(domain,xwidth,ywidth,depths):
    '''
    Builds a horizontally layered box like model. All layers are 
    tagged as 'intface_i' where i is the python style integer denoting
    that layer. For example, the free surface is tagged 'interface_0'.
    Volumes are similarly tagged as 'volume_i'.

    Parameters:
       domain :: output of Pycad.Design - It needs to be dim 3.\
       xwidth :: width in x direction in meters.\
       ywidth :: width in y direction in meters.\
       depth  :: depth to bottom of layer in meters.\

    One may save the domain using:
        # Output settings.\    
        domain.setScriptFileName(os.path.join(save_path,fname+".geo"))\
        domain.setMeshFileName(os.path.join(save_path,fname+".msh"))\
        findomain=fin.MakeDomain(domain) #  make the finley domain:\
    
        Create a file that can be read back in to python with ReadMesh.\
        findomain.write(os.path.join(save_path,fname+".fly"))\             
    '''
    
    #get number of layers
    if not hasattr(depths,"__len__"): depths = [ depths, ]
    ndepths=len(depths)

    if domain.getDim() != 3:
        raise TypeError("domain must be of dimension order 3.")        

    # Build the First Surface and add it to the domain
    fsuf,fsurl,fsurp=buildFreeSurface(xwidth,ywidth)
    domain.addItems(PropertySet('intface_%d'%(0),fsuf))        

    # Build each layer sequentially
    # Set up temporary variables.
    tsuf=fsuf; tsurl=fsurl; tsurp=fsurp
    # Loop through all depths.
    for i in range(0,ndepths):
        # Build new layer.
        tvol,tsuf,tsurl,tsurp=buildLayer(xwidth,ywidth,depths[i],\
                                     tsuf,tsurl,tsurp)
        # Add the new interface to the domain.
        domain.addItems(PropertySet('intface_%d'%(i+1),tsuf))        
        # Add the new volume/layer to the domain.
        domain.addItems(PropertySet('volume_%d'%i,tvol))
    
    return domain

