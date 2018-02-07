##############################################################################
#
# Copyright (c) 2009-2018 by The University of Queensland
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
from __future__ import division, print_function

__copyright__="""Copyright (c) 2009-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""


############################################################FILE HEADER
# example09n.py
# Create a simple 3D model for use in example09. This is the low res
# mesh for illustration purposes only.
#
#######################################################EXTERNAL MODULES
from esys.pycad import * #domain constructor
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.escript import mkDir, getMPISizeWorld
import os
import math
try:
    # This imports the rectangle domain function 
    from esys.finley import MakeDomain
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
        import sys
        print("This example will not run in an MPI world.")
        sys.exit(0)

if HAVE_FINLEY:
    # make sure path exists 
    save_path= os.path.join("data","example09c") 
    mkDir(save_path)

    ################################################ESTABLISHING PARAMETERS
    #Model Parameters
    origin=[0,0]                  #orign of model
    xwidth=300.0                      #width of model
    depth=-100.0                      #depth of model
    nintf=3                         #number of the interfaces
    lintf_depths=[-20,-40,-60]       #depth of interfaces
    rintf_depths=[-30,-50,-70]      #vertical displacement across fault
    fault_dip=40.0                  #dip of fault plane
    fault_atsurface=50.0             #location of fault at surface
    fault_width=10                 #apparent width of fault plane 

    element_size=1.

    ####################################################DOMAIN CONSTRUCTION
    x0=0.0+origin[0]
    y0=0.0+origin[1]
    #z=0.0+origin[2]
    fault_atsurface=fault_atsurface+origin[0]
    xwidth=xwidth+origin[0]
    depth=depth+origin[1]
    # Construction points, 4 vectors that descend from the surface with nintf+2 points.
    left_edge=[Point(x0,y0)]; 
    leftf_edge=[Point(fault_atsurface,y0)]; 
    rightf_edge=[Point(fault_atsurface+fault_width,y0)]; 
    right_edge=[Point(xwidth,y0)]; 

    for i in range(0,nintf):
        left_shift=(lintf_depths[i]-y0)/math.tan(fault_dip)
        right_shift=(rintf_depths[i]-y0)/math.tan(fault_dip)
        left_edge.append(Point(x0,lintf_depths[i]+origin[1]))
        leftf_edge.append(Point(fault_atsurface+left_shift,lintf_depths[i]+origin[1]))
        rightf_edge.append(Point(fault_atsurface+right_shift+fault_width,rintf_depths[i]+origin[1]))
        right_edge.append(Point(xwidth,rintf_depths[i]+origin[1]))

    left_edge.append(Point(x0,depth))
    leftf_edge.append(Point(fault_atsurface+(depth-y0)/math.tan(fault_dip),depth))
    rightf_edge.append(Point(fault_atsurface+fault_width+(depth-y0)/math.tan(fault_dip),depth))
    right_edge.append(Point(xwidth,depth))

    #Build lines
    lright=[]; nlright=[];
    lhright=[]; nlhright=[];
    lfright=[]; nlfright=[];
    lfhor=[]; nlfhor=[];
    lfleft=[]; nlfleft=[];
    lhleft=[]; nlhleft=[];
    lleft=[]; nlleft=[];

    #Build vertical lines
    for i in range(0,nintf+1):
        lleft.append(Line(left_edge[i],left_edge[i+1]))
        lfleft.append(Line(leftf_edge[i],leftf_edge[i+1]))
        lfright.append(Line(rightf_edge[i],rightf_edge[i+1]))
        lright.append(Line(right_edge[i],right_edge[i+1]))

    #Build horizontal lines
    for i in range(0,nintf+2):
        lhleft.append(Line(left_edge[i],leftf_edge[i]))
        lhright.append(Line(rightf_edge[i],right_edge[i]))
    lfhor.append(Line(leftf_edge[0],rightf_edge[0]))
    lfhor.append(Line(leftf_edge[nintf+1],rightf_edge[nintf+1]))

    #Build negative lines
    for i in range(nintf,-1,-1):
        nlleft.append(-lleft[i])
        nlfleft.append(-lfleft[i])
        nlfright.append(-lfright[i])
        nlright.append(-lright[i])
    for i in range(nintf+1,-1,-1):
        nlhleft.append(-lhleft[i])
        nlhright.append(-lhright[i])

    #Build curveloops
    lcurves=[]
    fcurves=[]
    rcurves=[]

    #Fault
    for i in range(0,nintf+1):
        fcurves.append(lfleft[i])
    fcurves.append(lfhor[1])
    for i in range(0,nintf+1):
        fcurves.append(nlfright[i])
    fcurves.append(-lfhor[0])
    fcurves=CurveLoop(*tuple(fcurves))

    #Left and Right Blocks
    for i in range(0,nintf+1):
        lcurves.append(CurveLoop(lleft[i],lhleft[i+1],nlfleft[nintf-i],nlhleft[nintf+1-i]))
        rcurves.append(CurveLoop(lfright[i],lhright[i+1],nlright[nintf-i],nlhright[nintf+1-i]))
        
    #Build Surfaces
    fsurf=PlaneSurface(fcurves)
    lsurf=[]
    rsurf=[]
    for i in range(0,nintf+1):
        lsurf.append(PlaneSurface(lcurves[i]))
        rsurf.append(PlaneSurface(rcurves[i]))

    d=Design(dim=2, element_size=element_size, order=2)

    d.addItems(PropertySet('fault',fsurf))
    for i in range(0,nintf+1):
        d.addItems(PropertySet('lblock%d'%i,lsurf[i]))
        d.addItems(PropertySet('rblock%d'%i,rsurf[i]))

    d.addItems(PropertySet('top',lhright[0],lfhor[0],lhleft[0],lhright[4],lhleft[4],lfhor[1]))

    d.setScriptFileName(os.path.join(save_path,"example09n.geo"))
    d.setMeshFileName(os.path.join(save_path,"example09n.msh"))
    #
    #  make the domain:
    #
    domain=MakeDomain(d)
    # mesh=ReadMesh(fileName) this is how to read the fly file into escript
    domain.write(os.path.join(save_path,"example09n.fly"))





