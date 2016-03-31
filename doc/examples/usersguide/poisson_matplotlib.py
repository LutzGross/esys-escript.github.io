##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import matplotlib
matplotlib.use('agg')    #For interactive use, you can comment out this line
#It's just here to make testing easier

import numpy
import pylab 

from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys.finley import Rectangle

#Testing whether we have a late enough version of matplotlib
try:
    matplotlib.mlab.griddata

    # TO keep the version distributed by openSuse happy
    interp='nn'
    try:
        from mpl_toolkits.natgrid import _natgrid
    except ImportError:
        interp='linear'
    # generate domain:
    mydomain = Rectangle(l0=1.,l1=1.,n0=40, n1=20)
    # define characteristic function of Gamma^D
    x = mydomain.getX()
    gammaD = whereZero(x[0])+whereZero(x[1])
    # define PDE and get its solution u
    mypde = Poisson(domain=mydomain)
    mypde.setValue(f=1,q=gammaD)
    u = mypde.getSolution()
    
    # interpolate u to a matplotlib grid:
    x_grid = numpy.linspace(0.,1.,50)
    y_grid = numpy.linspace(0.,1.,50)
    x=mydomain.getX()[0].toListOfTuples()
    y=mydomain.getX()[1].toListOfTuples()
    z=interpolate(u,mydomain.getX().getFunctionSpace()).toListOfTuples()
    z_grid = matplotlib.mlab.griddata(x,y,z,xi=x_grid,yi=y_grid,interp=interp )
    # interpolate u to a rectangular grid:
    matplotlib.pyplot.contourf(x_grid, y_grid, z_grid, 5)
    matplotlib.pyplot.savefig("u.png")
    # uncomment this line if you want to interact with a plot window
    #matplotlib.pyplot.show()

except AttributeError:
    print("Your version of matplotlib does not provide the griddata method.\nSkipping example.\n")
