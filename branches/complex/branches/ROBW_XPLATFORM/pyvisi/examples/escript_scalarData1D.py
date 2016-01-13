#!/usr/bin/env python


# $Id: escript_scalarData1D.py,v 1.3 2006/01/03 08:46:20 paultcochrane Exp $

from esys.escript import *
import esys.finley as fe

from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

online=True

vectorDomain = fe.Rectangle(9,10)
vectorFunctionSpace = escript.ContinuousFunction(vectorDomain)
domainData = vectorFunctionSpace.getX()

# plotting scalar data in a 1D array
scalarData1D = sin(domainData[0])

scene = Scene()
plot = LinePlot(scene)
plot.setData(scalarData1D)
if online: scene.render(pause=True)

scene.save(fname="escript_scalarData1D_linePlot.png", format="png")

# what other kinds of plot should I use for this kind of data??

# vim: expandtab shiftwidth=4:
