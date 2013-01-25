#!/usr/bin/env python


# $Id: escript_threeDVectorData3D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

brickDomain = bruce.Brick(9,9,9,10,10,10)
brickFunctionSpace=escript.ContinuousFunction(brickDomain)

vectorData3D = numarray.array([1.0,2.0,3.0])

# plotting 3D vectors in a 3D array
threeDVectorData3D = Data(vectorData3D, brickFunctionSpace, True)
scene = Scene()
plot = ArrowPlot3D(scene)
plot.setData(threeDVectorData3D)
if online: scene.render(pause=True)

scene.save(fname="escript_threeDVectorData3D_arrowPlot3D.png", format="png")

# vim: expandtab shiftwidth=4:
