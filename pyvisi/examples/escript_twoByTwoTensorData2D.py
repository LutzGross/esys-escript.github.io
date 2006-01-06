#!/usr/bin/env python


# $Id: escript_twoByTwoTensorData2D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from pyvisi import *
from pyvisi.renderers.vtk import *

tensorDomain = bruce.Rectangle(9,9,10,10)
tensorFunctionSpace=escript.ContinuousFunction(tensorDomain)

tensorData2x2 = numarray.array([[1.0,2.0],[3.0,4.0]])

# plotting 2x2 tensors in a 2D array
twoByTwoTensorData2D = Data(tensorData2x2, tensorFunctionSpace, True)
scene = Scene()
plot = EllipsoidPlot(scene)
plot.setData(twoByTwoTensorData2D)
scene.render(pause=True)

scene.save(fname="escript_twoByTwoTensorData2D_ellipsoidPlot.png", 
	format="png")

# vim: expandtab shiftwidth=4:
