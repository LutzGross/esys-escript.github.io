from __future__ import division, print_function
##############################################################################
#
# Copyright (c) 2009-2015 by University of Queensland
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

__copyright__="""Copyright (c) 2009-2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example05c.py
# Create either a 2D syncline or anticline model using pycad meshing 
# tools. The model steady state heat solution. How to create arrow or
# quiver plots.

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
from esys.escript.pdetools import Projector
from cblib import toRegGrid, subsample
import pylab as pl #Plotting package
import numpy as np

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

if HAVE_FINLEY:
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
    ###############################################################PLOTTING
    # show temperature:
    xi, yi, zi = toRegGrid(T, nx=50, ny=50)
    CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
    pl.clabel(CS, inline=1, fontsize=8)
    # show sub domains:
    tpg=np.array([p.getCoordinates() for p in tblockloop.getPolygon() ])
    pl.fill(tpg[:,0],tpg[:,1],'brown',label='2 W/m/k',zorder=-1000)
    bpg=np.array([p.getCoordinates() for p in bblockloop.getPolygon() ])
    pl.fill(bpg[:,0],bpg[:,1],'red',label='4 W/m/k',zorder=-1000)
    # show flux:
    xflux, flux=subsample(-kappa*grad(T), nx=20, ny=20)
    pl.quiver(xflux[:,0],xflux[:,1],flux[:,0],flux[:,1], angles='xy',color="white")
    # create plot
    pl.title("Heat Refraction across a clinal structure\n with heat flux.")
    pl.xlabel("Horizontal Displacement (m)")
    pl.ylabel("Depth (m)")
    pl.legend()
    pl.savefig(os.path.join(save_path,"flux.png"))
    print("Flux has been plotted  ...")
