# $Id: multiCurveLinePlot.py,v 1.1 2005/05/05 01:57:24 paultcochrane Exp $

"""
Example of plotting multiple curves with pyvisi 

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


# what plotting method are we using?
method = 'pyvisi'

# set up some data to plot
from Numeric import *

x = arange(0, 2*pi, 0.1, typecode=Float)
y1 = sin(x)
y2 = cos(x)
y3 = cos(x)**2

import plplot

# determine the min and max of x
xMin = min(x)
xMax = max(x)

# determine the global min and max of all the y's
yAll = concatenate( [y1, y2, y3] )

yMin = min(yAll)
yMax = max(yAll)

plplot.plsdev("xwin")
plplot.plinit()
plplot.plenv(xMin, xMax, yMin, yMax, 0, 1)
plplot.pllab("x", "y", "Example 2D plot")
plplot.plline(x, y1)
plplot.plline(x, y2)
plplot.plline(x, y3)
plplot.plend()

# to save as well, have to set everything up again, and replot
# save as png
plplot.plsdev("png")
plplot.plsfnam("multiCurveLinePlot.png")
plplot.plinit()
plplot.plenv(xMin, xMax, yMin, yMax, 0, 1)
plplot.pllab("x", "y", "Example 2D plot")
plplot.plline(x, y1)
plplot.plline(x, y2)
plplot.plline(x, y3)
plplot.plend()

