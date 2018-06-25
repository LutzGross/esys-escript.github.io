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
# example06.py
# Create a 3 block fault and overburden style model in 2d using pycad
# meshing tools.

#######################################################EXTERNAL MODULES
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
from esys.pycad import *
from esys.pycad.gmsh import Design
from esys.escript import *
import numpy as np
import pylab as pl #Plotting package
from cblib import toRegGrid, subsample, HAVE_NATGRID
from esys.escript.unitsSI import *
from esys.escript.linearPDEs import LinearPDE
import os, sys

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
    # where to put output files
    save_path= os.path.join("data","example06")
    mkDir(save_path)

    Ttop=20*Celsius # temperature at the top
    qin=300.*Milli*W/(m*m) #our heat source temperature is now zero

    ################################################ESTABLISHING PARAMETERS
    # Overall Domain
    p0=Point(0.0,        0.0, 0.0)
    p1=Point(0.0,    -6000.0, 0.0)
    p2=Point(5000.0, -6000.0, 0.0)
    p3=Point(5000.0,     0.0, 0.0)

    ####################################################DOMAIN CONSTRUCTION
    l01=Line(p0, p1)
    l12=Line(p1, p2)
    l23=Line(p2, p3)
    l30=Line(p3, p0)

    # Generate Material Boundary
    p4=Point(0.0,    -2400.0, 0.0)
    p5=Point(2000.0, -2400.0, 0.0)
    p6=Point(3000.0, -6000.0, 0.0)
    p7=Point(5000.0, -2400.0, 0.0)

    # Create TOP BLOCK
    tbl1=Line(p0,p4)
    tbl2=Line(p4,p5)
    tbl3=Line(p5,p7)
    tbl4=Line(p7,p3)
    tblockloop = CurveLoop(tbl1,tbl2,tbl3,tbl4,l30)
    tblock = PlaneSurface(tblockloop)

    # Create BOTTOM BLOCK LEFT
    bbll1=Line(p4,p1)
    bbll2=Line(p1,p6)
    bbll3=Line(p6,p5)
    bbll4=-tbl2
    bblockloopl = CurveLoop(bbll1,bbll2,bbll3,bbll4)
    bblockl = PlaneSurface(bblockloopl)

    # Create BOTTOM BLOCK RIGHT
    bbrl1=Line(p6,p2)
    bbrl2=Line(p2,p7)
    bbrl3=-tbl3
    bbrl4=-bbll3
    bblockloopr = CurveLoop(bbrl1,bbrl2,bbrl3,bbrl4)
    bblockr = PlaneSurface(bblockloopr)

    #############################################EXPORTING MESH FOR ESCRIPT
    # Create a Design which can make the mesh
    d=Design(dim=2, element_size=200)
    # Add the subdomains and flux boundaries.
    d.addItems(PropertySet("top",tblock),\
                           PropertySet("bottomleft",bblockl),\
                           PropertySet("bottomright",bblockr),\
                           PropertySet("linebottom",bbll2, bbrl1))

    # Create the geometry, mesh and Escript domain
    d.setScriptFileName(os.path.join(save_path,"example06.geo"))

    d.setMeshFileName(os.path.join(save_path,"example06.msh"))
    domain=MakeDomain(d)
    print("Domain has been generated ...")

    # set up kappa (thermal conductivity across domain) using tags
    kappa=Scalar(0,Function(domain))
    kappa.setTaggedValue("top",2.0)
    kappa.setTaggedValue("bottomleft",10.0)
    kappa.setTaggedValue("bottomright",6.0)
    ##############################################################SOLVE PDE
    mypde=LinearPDE(domain)
    mypde.getSolverOptions().setVerbosityOn()
    mypde.setSymmetryOn()
    mypde.setValue(A=kappa*kronecker(domain))
    x=Solution(domain).getX()
    mypde.setValue(q=whereZero(x[1]-sup(x[1])),r=Ttop)
    qS=Scalar(0,FunctionOnBoundary(domain))
    qS.setTaggedValue("linebottom",qin)
    mypde.setValue(y=qS)
    print("PDE has been generated ...")
    ###########################################################GET SOLUTION
    T=mypde.getSolution()
    print("PDE has been solved ...")
    ###############################################################PLOTTING
    # show temperature:
    xi, yi, zi = toRegGrid(T, nx=50, ny=50)
    CS = pl.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
    pl.clabel(CS, inline=1, fontsize=8)
    # show sub domains:
    tpg=np.array([p.getCoordinates() for p in tblockloop.getPolygon() ])
    pl.fill(tpg[:,0],tpg[:,1],'brown',label='2 W/m/k',zorder=-1000)
    bpgr=np.array([p.getCoordinates() for p in bblockloopr.getPolygon() ])
    pl.fill(bpgr[:,0],bpgr[:,1],'orange',label='6 W/m/k',zorder=-1000)
    bpgl=np.array([p.getCoordinates() for p in bblockloopl.getPolygon() ])
    pl.fill(bpgl[:,0],bpgl[:,1],'red',label='10 W/m/k',zorder=-1000)
    # show flux:
    xflux, flux=subsample(-kappa*grad(T), nx=20, ny=20)
    pl.quiver(xflux[:,0],xflux[:,1],flux[:,0],flux[:,1], angles='xy',color="white")
    # create plot
    pl.title("Heat Refraction across an anisotropic structure\n with flux vectors.")
    pl.xlabel("Horizontal Displacement (m)")
    pl.ylabel("Depth (m)")
    pl.legend()
    pl.savefig(os.path.join(save_path,"flux.png"))
    print("Flux has been plotted  ...")
