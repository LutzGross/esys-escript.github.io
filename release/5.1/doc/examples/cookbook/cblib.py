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
Additional routines using matplotlib for cookbook examples.
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""
from esys.escript import inf,sup
from esys.escript.pdetools import Locator
import numpy as np
import pylab as pl

try:
    from mpl_toolkits.natgrid import _natgrid
    HAVE_NATGRID=True
except ImportError:
    HAVE_NATGRID=False

def toXYTuple(coords):
    """
    extracts the X and Y coordinates as two ```numpy`` arrays from the escript coordinates ```coords``` as produced by a ``.getX`` call.
    """
    coords = np.array(coords.toListOfTuples()) #convert to array
    coordX = coords[:,0]; coordY = coords[:,1] #X and Y components.
    return coordX,coordY

def toRegGrid(u, nx=50, ny=50):
   """
   returns a nx x ny grid representation of the escript object u
   """
   xx=u.getDomain().getX()     
   x=u.getFunctionSpace().getX()     
   coordX, coordY = toXYTuple(x)
   utemp = u.toListOfTuples()
   # create regular grid
   xi = np.linspace(inf(xx[0]),sup(xx[0]),nx)
   yi = np.linspace(inf(xx[1]),sup(xx[1]),ny)

   # interpolate u to grid
   zi = pl.matplotlib.mlab.griddata(coordX,coordY,utemp,xi,yi, interp='linear')
   return xi, yi, zi

def subsample(u, nx=50, ny=50):
    """
    subsamples ```u``` over an ```nx``` x ```ny``` grid
    and returns ``numpy`` arrays for the values and locations
    used for subsampling.
    """
    xx=u.getDomain().getX()  # points of the domain
    x0=inf(xx[0])         
    y0=inf(xx[1])
    dx = (sup(xx[0])-x0)/nx # x spacing
    dy = (sup(xx[1])-y0)/ny # y spacing
    grid = [ ] 
    for j in range(0,ny-1):
        for i in range(0,nx-1):
               grid.append([x0+dx/2+dx*i,y0+dy/2+dy*j])
    uLoc = Locator(u.getFunctionSpace(),grid)
    subu= uLoc(u) # get data of u at sample points closests to grid points
    usublocs = uLoc.getX() #returns actual locations from data
    return np.array(usublocs), np.array(subu)
    
