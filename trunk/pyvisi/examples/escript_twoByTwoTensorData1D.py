#!/usr/bin/env python


# $Id: escript_twoByTwoTensorData1D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from pyvisi import *
from pyvisi.renderers.vtk import *

vectorDomain = bruce.Rectangle(9,1,10,1)
vectorFunctionSpace=escript.ContinuousFunction(vectorDomain)

tensorData2x2 = numarray.array([[1.0,2.0],[3.0,4.0]])

# plotting 2x2 tensors in a 1D array
twoByTwoTensorData1D = Data(tensorData2x2, vectorFunctionSpace, True)
scene = Scene()
plot = EllipsoidPlot(scene)
plot.setData(twoByTwoTensorData1D)
scene.render(pause=True)

scene.save(fname="escript_twoByTwoTensorData1D_ellipsoidPlot.png", 
	format="png")

# vim: expandtab shiftwidth=4:
