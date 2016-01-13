#!/usr/bin/env python


# $Id: escript_scalarData2D.py,v 1.3 2006/01/03 08:46:20 paultcochrane Exp $

from esys.escript import *
import esys.finley as fe

import numarray

# from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

online=True

tensorDomain = fe.Rectangle(9,9)
tensorFunctionSpace=escript.ContinuousFunction(tensorDomain)
domainData = tensorFunctionSpace.getX()

# plotting scalar data in a 2D array
scalarData2D = sin(domainData[0])

scene = Scene()
plot = ContourPlot(scene)
plot.setData(scalarData2D)
if online: scene.render(pause=True)

scene.save(fname="escript_scalarData2D_contourPlot.png", format="png")

scene = Scene()
plot = SurfacePlot(scene)
plot.setData(-scalarData2D)
if online: scene.render(pause=True)
plot.setData(-scalarData2D)
if online: scene.render(pause=True)

# scene.save(fname="escript_scalarData2D_surfacePlot.png", format="png")

# add SurfacePlot, MeshPlot etc here

# vim: expandtab shiftwidth=4:
