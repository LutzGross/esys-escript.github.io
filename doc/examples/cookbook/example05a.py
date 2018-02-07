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
# example05a.py
# Create either a 2D syncline or anticline model using pycad meshing 
# tools. Then model steady state heat solution.

#######################################################EXTERNAL MODULES
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.pycad import * #domain constructor
from esys.pycad.gmsh import Design #Finite Element meshing package
import os #file path tool
from math import * # math package
from esys.escript import *
from esys.escript.unitsSI import *
from esys.escript.linearPDEs import LinearPDE
from cblib import toRegGrid, HAVE_NATGRID
import pylab as pl #Plotting package

try:
    # This imports the rectangle domain function 
    from esys.finley import MakeDomain#Converter for escript
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
        import sys
        print("This example will not run in an MPI world.")
        sys.exit(0)

if not HAVE_NATGRID:
    print("This example requires that natgrid is available to matplotlib")

if HAVE_FINLEY and HAVE_NATGRID:
    #################################################ESTABLISHING VARIABLES
    #set modal to 1 for a syncline or -1 for an anticline structural 
    #configuration
    modal=-1

    # the folder to put our outputs in, leave blank "" for script path - 
    # note this folder path must exist to work
    save_path= os.path.join("data","example05") 
    mkDir(save_path)

    ################################################ESTABLISHING PARAMETERS
    #Model Parameters
    width=5000.0*m   #width of model
    depth=-6000.0*m  #depth of model
    Ttop=20*K       # top temperature
    qin=70*Milli*W/(m*m) # bottom heat influx

    sspl=51 #number of discrete points in spline
    dsp=width/(sspl-1) #dx of spline steps for width
    dep_sp=2500.0*m #avg depth of spline
    h_sp=1500.0*m #heigh of spline
    orit=-1.0 #orientation of spline 1.0=>up -1.0=>down

    ####################################################DOMAIN CONSTRUCTION
    # Domain Corners
    p0=Point(0.0,      0.0, 0.0)
    p1=Point(0.0,    depth, 0.0)
    p2=Point(width, depth, 0.0)
    p3=Point(width,   0.0, 0.0)

    # Generate Material Boundary
    x=[ Point(i*dsp\
        ,-dep_sp+modal*orit*h_sp*cos(pi*i*dsp/dep_sp+pi))\
         for i in range(0,sspl)\
      ]
    mysp = Spline(*tuple(x))
    # Start and end of material boundary.
    x1=mysp.getStartPoint()
    x2=mysp.getEndPoint()

    #  Create TOP BLOCK
    # lines
    tbl1=Line(p0,x1)
    tbl2=mysp
    tbl3=Line(x2,p3)
    l30=Line(p3, p0)
    # curve
    tblockloop = CurveLoop(tbl1,tbl2,tbl3,l30)
    # surface
    tblock = PlaneSurface(tblockloop)


    # Create BOTTOM BLOCK
    # lines
    bbl1=Line(x1,p1)
    bbl3=Line(p2,x2)
    bbl4=-mysp
    l12=Line(p1, p2)
    # curve
    bblockloop = CurveLoop(bbl1,l12,bbl3,bbl4)
    # surface
    bblock = PlaneSurface(bblockloop)

    ################################################CREATE MESH FOR ESCRIPT
    # Create a Design which can make the mesh
    d=Design(dim=2, element_size=200)
    # Add the subdomains and flux boundaries.
    d.addItems(PropertySet("top",tblock),PropertySet("bottom",bblock),\
                                         PropertySet("linebottom",l12))
    # Create the geometry, mesh and Escript domain
    d.setScriptFileName(os.path.join(save_path,"example05.geo"))
    d.setMeshFileName(os.path.join(save_path,"example05.msh"))
    domain=MakeDomain(d, optimizeLabeling=True)
    print("Domain has been generated ...")
    ##############################################################SOLVE PDE
    mypde=LinearPDE(domain)
    mypde.getSolverOptions().setVerbosityOn()
    mypde.setSymmetryOn()
    kappa=Scalar(0,Function(domain))
    kappa.setTaggedValue("top",2.0*W/m/K)
    kappa.setTaggedValue("bottom",4.0*W/m/K)
    mypde.setValue(A=kappa*kronecker(domain))
    x=Solution(domain).getX()
    mypde.setValue(q=whereZero(x[1]-sup(x[1])),r=Ttop)
    qS=Scalar(0,FunctionOnBoundary(domain))
    qS.setTaggedValue("linebottom",qin)
    mypde.setValue(y=qS)
    print("PDE has been generated ...")
    ###########################################################GET SOLUTION
    T=mypde.getSolution()
    print("PDE has been solved  ...")

    #######################################################################
    xi, yi, zi = toRegGrid(T, nx=50, ny=50)
    pl.matplotlib.pyplot.autumn()
    pl.contourf(xi,yi,zi,10)
    pl.xlabel("Horizontal Displacement (m)")
    pl.ylabel("Depth (m)")
    pl.savefig(os.path.join(save_path,"Tcontour.png"))
    print("Solution has been plotted  ...")
