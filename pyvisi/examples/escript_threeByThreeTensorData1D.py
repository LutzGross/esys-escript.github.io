#!/usr/bin/env python


# $Id: escript_threeByThreeTensorData1D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from pyvisi import *
from pyvisi.renderers.vtk import *

vectorDomain = bruce.Rectangle(9,1,10,1)
vectorFunctionSpace=escript.ContinuousFunction(vectorDomain)

tensorData3x3 = numarray.array([[1.0,2.0,3.0],[3.0,4.0,5.0],[5.0,6.0,7.0]])

# plotting 3x3 tensors in a 1D array
threeByThreeTensorData1D = Data(tensorData3x3, vectorFunctionSpace, True)
scene = Scene()
plot = EllipsoidPlot(scene)
plot.setData(threeByThreeTensorData1D)
scene.render(pause=True)

scene.save(fname="escript_threeByThreeTensorData1D_ellipsoidPlot.png",
	format="png")

# vim: expandtab shiftwidth=4:
