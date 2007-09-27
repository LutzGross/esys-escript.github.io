#!/usr/bin/env python

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
