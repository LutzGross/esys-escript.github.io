#!/usr/bin/env python


# $Id: escript_threeDVectorData1D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

online=False

vectorDomain = bruce.Rectangle(9,1,10,1)
vectorFunctionSpace=escript.ContinuousFunction(vectorDomain)

vectorData3D = numarray.array([1.0,2.0,3.0])

# plotting 3D vectors in a 1D array
threeDVectorData1D = Data(vectorData3D, vectorFunctionSpace, True)
scene = Scene()
plot = ArrowPlot3D(scene)
plot.setData(threeDVectorData1D)
if online: scene.render(pause=True)

scene.save(fname="escript_threeDVectorData1D_arrowPlot3D.png", format="png")

# vim: expandtab shiftwidth=4:
