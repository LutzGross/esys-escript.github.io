
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys.finley import Rectangle
import numpy
import matplotlib

matplotlib.use('agg')	#For interactive use, you can comment out this line
#It's just here to make testing easier

import pylab 

#Testing whether we have a late enough version of matplotlib
try:
	matplotlib.mlab.griddata
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
	z_grid = matplotlib.mlab.griddata(x,y,z,xi=x_grid,yi=y_grid )
	# interpolate u to a rectangular grid:
	matplotlib.pyplot.contourf(x_grid, y_grid, z_grid, 5)
	matplotlib.pyplot.savefig("u.png")
	# uncommend this line if you want to interact with a plot window
	#matplotlib.pyplot.show()

except AttributeError:
	print "Your version of matplotlib does not provide the griddata method.\nSkipping example.\n"
