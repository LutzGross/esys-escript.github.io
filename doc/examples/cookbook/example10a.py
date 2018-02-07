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
# example10a.py
# Model of gravitational Potential for a gravity POLE.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
import matplotlib
matplotlib.use('agg') #It's just here for automated testing

from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.weipa import saveVTK # This imports the VTK file saver from weipa
import os, sys #This package is necessary to handle saving our data.
from math import pi, sqrt, sin, cos

from esys.escript.pdetools import Projector

from cblib import toRegGrid, HAVE_NATGRID
import pylab as pl #Plotting package
import numpy as np

try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
    print("This example will not run in an MPI world.")
    sys.exit(0)

if not HAVE_NATGRID:
    print("This example requires that natgrid is available to matplotlib")

if HAVE_FINLEY and HAVE_NATGRID:
    #################################################ESTABLISHING VARIABLES
    #Domain related.
    mx = 5000*m #meters - model length
    my = -5000*m #meters - model width
    ndx = 100 # mesh steps in x direction 
    ndy = 100 # mesh steps in y direction - one dimension means one element
    #PDE related
    rho=200.0
    rholoc=[2500,-2500]
    G=6.67300*10E-11

    ################################################ESTABLISHING PARAMETERS
    #the folder to put our outputs in, leave blank "" for script path 
    save_path= os.path.join("data","example10")
    #ensure the dir exists
    mkDir(save_path)

    ####################################################DOMAIN CONSTRUCTION
    domain = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
    x=Solution(domain).getX()
    mask=wherePositive(10-length(x-rholoc))
    rho=rho*mask
    kro=kronecker(domain)

    q=whereZero(x[1]-my)+whereZero(x[1])+whereZero(x[0])+whereZero(x[0]-mx)
    ###############################################ESCRIPT PDE CONSTRUCTION

    mypde=LinearPDE(domain)
    mypde.setValue(A=kro,Y=4.*3.1415*G*rho)
    mypde.setValue(q=q,r=0)
    mypde.setSymmetryOn()
    sol=mypde.getSolution()

    g_field=grad(sol) #The gravitational acceleration g.
    g_fieldz=g_field*[0,1] #The vertical component of the g field.
    gz=length(g_fieldz) #The magnitude of the vertical component.
    # Save the output to file.
    saveVTK(os.path.join(save_path,"ex10a.vtu"),\
            grav_pot=sol,g_field=g_field,g_fieldz=g_fieldz,gz=gz)

    ##################################################REGRIDDING & PLOTTING


    xi, yi, zi = toRegGrid(sol, nx=50, ny=50)
    pl.matplotlib.pyplot.autumn()
    pl.contourf(xi,yi,zi,10)
    pl.xlabel("Horizontal Displacement (m)")
    pl.ylabel("Depth (m)")
    pl.savefig(os.path.join(save_path,"Ucontour.png"))
    print("Solution has been plotted  ...")

    cut=int(len(xi)//2)

    pl.clf()

    r=np.linspace(0.0000001,mx/2,100)   # starting point would be 0 but that would cause division by zero later
    m=2*pl.pi*10*10*200*-G/(r*r)

    pl.plot(xi,zi[:,cut])
    #pl.plot(r+2500,m)
    pl.title("Potential Profile")
    pl.xlabel("Horizontal Displacement (m)")
    pl.ylabel("Potential")
    pl.savefig(os.path.join(save_path,"Upot00.png"))

    out=np.array([xi,zi[:,cut]])
    pl.savetxt('profile1.asc',out.transpose())
    pl.clf()
