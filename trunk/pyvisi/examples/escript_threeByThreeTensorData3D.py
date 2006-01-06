#!/usr/bin/env python


# $Id: escript_threeByThreeTensorData3D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from pyvisi import *
from pyvisi.renderers.vtk import *

brickDomain = bruce.Brick(9,9,9,10,10,10)
brickFunctionSpace=escript.ContinuousFunction(brickDomain)

tensorData3x3 = numarray.array([[1.0,2.0,3.0],[3.0,4.0,5.0],[5.0,6.0,7.0]])

# plotting 3x3 tensors in a 3D array
threeByThreeTensorData3D = Data(tensorData3x3, brickFunctionSpace, True)
scene = Scene()
plot = EllipsoidPlot(scene)
plot.setData(threeByThreeTensorData3D)
scene.render(pause=True)

scene.save(fname="escript_threeByThreeTensorData3D_ellipsoidPlot.png", 
	format="png")

# vim: expandtab shiftwidth=4:
