#!/usr/bin/env python


# $Id: escript_scalarData3D.py,v 1.3 2006/01/03 08:46:20 paultcochrane Exp $

from esys.escript import *
from esys import bruce

import numarray

from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

online=False

brickDomain = bruce.Brick(9,9,9,10,10,10)
brickFunctionSpace = escript.ContinuousFunction(brickDomain)

domainData = brickFunctionSpace.getX()

# plotting scalar data in a 3D array
scalarData3D = sin(domainData[0])

scene = Scene()
plot = IsosurfacePlot(scene)
plot.setData(scalarData3D)
if online: scene.render(pause=True)

scene.save(fname="escript_scalarData3D_isosurfacePlot.png", format="png")

# vim: expandtab shiftwidth=4:
