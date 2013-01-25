#!/usr/bin/env python


# $Id: escript_threeByThreeTensorData3D.py,v 1.3 2006/01/05 00:13:50 paultcochrane Exp $
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
from esys import bruce

import numarray

from esys.pyvisi import *
from esys.pyvisi.renderers.vtk import *

online=False

brickDomain = bruce.Brick(9,9,9,10,10,10)
brickFunctionSpace=escript.ContinuousFunction(brickDomain)

tensorData3x3 = numarray.array([[1.0,2.0,3.0],[3.0,4.0,5.0],[5.0,6.0,7.0]])

# plotting 3x3 tensors in a 3D array
threeByThreeTensorData3D = Data(tensorData3x3, brickFunctionSpace, True)
scene = Scene()
plot = EllipsoidPlot(scene)
plot.setData(threeByThreeTensorData3D)
if online: scene.render(pause=True)

scene.save(fname="escript_threeByThreeTensorData3D_ellipsoidPlot.png", 
	format="png")

# vim: expandtab shiftwidth=4:
