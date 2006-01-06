#!/usr/bin/env python


# $Id: escript_threeByThreeTensorData2D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from pyvisi import *
from pyvisi.renderers.vtk import *

tensorDomain = bruce.Rectangle(9,9,10,10)
tensorFunctionSpace=escript.ContinuousFunction(tensorDomain)

tensorData3x3 = numarray.array([[1.0,2.0,3.0],[3.0,4.0,5.0],[5.0,6.0,7.0]])

# plotting 3x3 tensors in a 2D array
threeByThreeTensorData2D = Data(tensorData3x3, tensorFunctionSpace, True)
scene = Scene()
plot = EllipsoidPlot(scene)
plot.setData(threeByThreeTensorData2D)
scene.render(pause=True)

scene.save(fname="escript_threeByThreeTensorData2D_ellipsoidPlot.png",
	format="png")

# vim: expandtab shiftwidth=4:
