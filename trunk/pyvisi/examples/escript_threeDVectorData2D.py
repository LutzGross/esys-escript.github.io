#!/usr/bin/env python


# $Id: escript_threeDVectorData2D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from pyvisi import *
from pyvisi.renderers.vtk import *

tensorDomain = bruce.Rectangle(9,9,10,10)
tensorFunctionSpace=escript.ContinuousFunction(tensorDomain)

vectorData3D = numarray.array([1.0,2.0,3.0])

# plotting 3D vectors in a 2D array
threeDVectorData2D = Data(vectorData3D, tensorFunctionSpace, True)
scene = Scene()
plot = ArrowPlot3D(scene)
plot.setData(threeDVectorData2D)
scene.render(pause=True)

scene.save(fname="escript_threeDVectorData2D_arrowPlot3D.png", format="png")

# vim: expandtab shiftwidth=4:
