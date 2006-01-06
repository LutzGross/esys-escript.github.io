#!/usr/bin/env python


# $Id: escript_scalarData1D.py,v 1.3 2006/01/03 08:46:20 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from pyvisi import *
from pyvisi.renderers.vtk import *

vectorDomain = bruce.Rectangle(9,1,10,1)
vectorFunctionSpace = escript.ContinuousFunction(vectorDomain)
domainData = vectorFunctionSpace.getX()

# plotting scalar data in a 1D array
scalarData1D = sin(domainData[0])

scene = Scene()
plot = LinePlot(scene)
plot.setData(scalarData1D)
scene.render(pause=True)

scene.save(fname="escript_scalarData1D_linePlot.png", format="png")

# what other kinds of plot should I use for this kind of data??

# vim: expandtab shiftwidth=4:
