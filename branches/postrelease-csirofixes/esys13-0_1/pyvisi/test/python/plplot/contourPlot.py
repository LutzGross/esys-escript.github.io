# $Id: contourPlot.py,v 1.1 2005/11/08 05:54:01 paultcochrane Exp $

"""
Example of contour plotting with pyvisi 

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


# set up some data to plot
from Numeric import *

# the x and y axes
x = arange(-2,2,0.2, typecode=Float)
y = arange(-2,3,0.2, typecode=Float)

# pick some interesting function to generate the data in the third dimension
# this is the one used in the matlab docs: z = x*exp(-x^2-y^2)
z = zeros((len(x),len(y)), typecode=Float)

# boy do *I* feel old fashioned writing it this way
# surely there's another way to do it: - something to do later
for i in range(len(x)):
    for j in range(len(y)):
	z[i,j] = x[i]*exp(-x[i]*x[i] - y[j]*y[j])

import plplot

# determine the min and max of x
xMin = min(x)
xMax = max(x)

yMin = min(y)
yMax = max(y)

plplot.plsdev("xwin")
plplot.plinit()
plplot.plenv(xMin, xMax, yMin, yMax, 0, 1)
plplot.pllab("x", "y", "Example shaded contour plot")
plshades(zz, shedge, fill_width, 1, pltr1, xg1, yg1)

zmin = min(zz.flat)
zmax = max(zz.flat)

clevel = zmin + (zmax - zmin) * (arrayrange(NS)+0.5)/NS
shedge = zmin + (zmax - zmin) * (arrayrange(NS+1))/NS

plplot.plend()

# to save as well, have to set everything up again, and replot
# save as png
plplot.plsdev("png")
plplot.plsfnam("contourPlot.png")
plplot.plinit()
plplot.plenv(xMin, xMax, yMin, yMax, 0, 1)
plplot.pllab("x", "y", "Example shaded contour plot")
plplot.plline(x, y1)
plplot.plend()

# vim: expandtab shiftwidth=4:
