##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import matplotlib
matplotlib.use('agg')    #For interactive use, you can comment out this line
#It's just here to make testing easier

import numpy
import pylab
try:
    import scipy.interpolate
    HAVE_SCIPY=True
except:
    HAVE_SCIPY=False

from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys.finley import Rectangle

#Testing whether we have a late enough version of matplotlib
if HAVE_SCIPY:
    try:
        interp = 'linear'
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
        z_grid = scipy.interpolate.griddata((x,y),z,(x_grid[None,:],y_grid[:,None]),interp)
        # interpolate u to a rectangular grid:
        matplotlib.pyplot.contourf(x_grid, y_grid, z_grid, 5)
        matplotlib.pyplot.savefig("u.png")
        # uncomment this line if you want to interact with a plot window
        #matplotlib.pyplot.show()

    except AttributeError:
        print("Your version of matplotlib does not provide the griddata method.\nSkipping example.\n")
else:
    print("This example requires scipy")
