##############################################################################
#
# Copyright (c) 2009-2017 by The University of Queensland
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
__copyright__="""Copyright (c) 2009-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example04b.py
# Create either a 2D mesh for a rectangle model using pycad meshing
# tools.
#
#######################################################EXTERNAL MODULES
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.pycad import * #domain constructor
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.escript import *
from esys.escript.unitsSI import *
from esys.escript.linearPDEs import LinearPDE
import pylab as pl #Plotting package
from cblib import toRegGrid, HAVE_NATGRID
import os
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

try:
    from mpl_toolkits.natgrid import _natgrid
    HAVE_NATGRID=True
except ImportError:
    HAVE_NATGRID=False

if not HAVE_NATGRID:
    print("This example requires that natgrid is available to matplotlib")

if HAVE_FINLEY and HAVE_NATGRID:
    # make sure path exists
    save_path= os.path.join("data","example04")
    mkDir(save_path)

    ################################################ESTABLISHING PARAMETERS
    #Model Parameters
    width=5000.0*m   #width of model
    depth=-6000.0*m  #depth of model
    kappa=2.0*W/m/K   # watts/m.Kthermal conductivity  
    Ttop=20*K       # top temperature
    qin=70*Milli*W/(m*m) # bottom heat influx

    ####################################################DOMAIN CONSTRUCTION
    # Domain Corners
    p0=Point(0.0,      0.0, 0.0)
    p1=Point(0.0,    depth, 0.0)
    p2=Point(width, depth, 0.0)
    p3=Point(width,   0.0, 0.0)
    # Join corners in anti-clockwise manner.
    l01=Line(p0, p1)
    l12=Line(p1, p2)
    l23=Line(p2, p3)
    l30=Line(p3, p0)
    # Join line segments to create domain boundary.
    c=CurveLoop(l01, l12, l23, l30)
    # surface
    rec = PlaneSurface(c)

    #############################################EXPORTING MESH FOR ESCRIPT
    # Create a Design which can make the mesh
    d=Design(dim=2, element_size=200*m)
    # Add the subdomains and flux boundaries.
    d.addItems(rec, PropertySet("linebottom",l12))
    d.addItems(l01, l23, l30) # just in case we need them
    #############################################MAKE THE DOMAIN
    domain=MakeDomain(d, optimizeLabeling=True)
    print("Domain has been generated ...")
    ##############################################################SOLVE PDE
    mypde=LinearPDE(domain)
    mypde.getSolverOptions().setVerbosityOn()
    mypde.setSymmetryOn()
    mypde.setValue(A=kappa*kronecker(domain))
    x=Solution(domain).getX()
    mypde.setValue(q=whereZero(x[1]-sup(x[1])),r=Ttop)
    qS=Scalar(0,FunctionOnBoundary(domain))
    qS.setTaggedValue("linebottom",qin)
    mypde.setValue(y=-qS)
    print("PDE has been generated ...")
    ###########################################################GET SOLUTION
    T=mypde.getSolution()
    print("PDE has been solved  ...")
    ###########################################################
    xi, yi, zi = toRegGrid(T, nx=50, ny=50)
    pl.matplotlib.pyplot.autumn()
    pl.contourf(xi,yi,zi,10)
    pl.xlabel("Horizontal Displacement (m)")
    pl.ylabel("Depth (m)")
    pl.savefig(os.path.join(save_path,"example04.png"))
    print("Solution has been plotted  ...")
