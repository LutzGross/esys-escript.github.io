#!/usr/bin/env python


# $Id: escript_scalarData2D.py,v 1.3 2006/01/03 08:46:20 paultcochrane Exp $
"""
@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"

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
