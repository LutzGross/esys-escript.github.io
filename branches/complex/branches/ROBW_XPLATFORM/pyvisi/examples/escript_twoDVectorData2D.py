#!/usr/bin/env python


# $Id: escript_twoDVectorData2D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

tensorDomain = bruce.Rectangle(9,9,10,10)
tensorFunctionSpace=escript.ContinuousFunction(tensorDomain)

vectorData2D = numarray.array([1.0,2.0])

# plotting 2D vectors in a 2D array
twoDVectorData2D = Data(vectorData2D, tensorFunctionSpace, True)
scene = Scene()
plot = ArrowPlot(scene)
plot.setData(twoDVectorData2D)
if online: scene.render(pause=True)

scene.save(fname="escript_twoDVectorData2D_arrowPlot.png", format="png")

# vim: expandtab shiftwidth=4:
