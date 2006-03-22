#!/usr/bin/env python


# $Id: escript_twoByTwoTensorData3D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

brickDomain = bruce.Brick(9,9,9,10,10,10)
brickFunctionSpace=escript.ContinuousFunction(brickDomain)

tensorData2x2 = numarray.array([[1.0,2.0],[3.0,4.0]])

# plotting 2x2 tensors in a 3D array
twoByTwoTensorData3D = Data(tensorData2x2, brickFunctionSpace, True)
scene = Scene()
plot = EllipsoidPlot(scene)
plot.setData(twoByTwoTensorData3D)
if online: scene.render(pause=True)

scene.save(fname="escript_twoByTwoTensorData3D_ellipsoidPlot.png", 
	format="png")

# vim: expandtab shiftwidth=4:
