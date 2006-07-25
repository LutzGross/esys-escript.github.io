# $Id: surfacePlot.py,v 1.1 2005/11/08 07:28:25 paultcochrane Exp $

"""
Example of plotting surfaces

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

# determine the min and max of x, y, and z in world coordinates
xMin = min(x)
xMax = max(x)

yMin = min(y)
yMax = max(y)

zMin = min(z.flat)
zMax = max(z.flat)

# min and max of x and y variables in normalised coordinates
# (these are values recommended by plplot in an example)
xMin2D = -2.5
xMax2D = 2.5

yMin2D = -2.5
yMax2D = 4.0

# sides of box in normalised coordinates
# (these are values recommended by plplot in an example)
basex = 2.0
basey = 4.0
height = 3.0

# angle to view box
alt = 45.0
az = 30.0

side = 1
opt = 3  # plots a net of lines

plplot.plsdev("xwin")
plplot.plinit()
plplot.plenv(xMin2D, xMax2D, yMin2D, yMax2D, 0, -2)
plplot.plw3d(basex, basey, height, 
	xMin, xMax, yMin, yMax, zMin, zMax, 
	alt, az)
plplot.plmtex("t", 1.0, 0.5, 0.5, "Example surface plot")
plplot.plbox3("bnstu", "x axis", 0.0, 0, 
	"bnstu", "y axis", 0.0, 0, 
	"bcdmnstuv", "z axis", 0.0, 0)
plplot.plsurf3d(x, y, z, 0, ())
plplot.plend()

# to save as well, have to set everything up again, and replot
# save as png
plplot.plsdev("png")
plplot.plsfnam("surfacePlot.png")
plplot.plinit()
plplot.plenv(xMin2D, xMax2D, yMin2D, yMax2D, 0, -2)
plplot.plw3d(basex, basey, height, 
	xMin, xMax, yMin, yMax, zMin, zMax, 
	alt, az)
plplot.plmtex("t", 1.0, 0.5, 0.5, "Example surface plot")
plplot.plbox3("bnstu", "x axis", 0.0, 0, 
	"bnstu", "y axis", 0.0, 0, 
	"bcdmnstuv", "z axis", 0.0, 0)
plplot.plsurf3d(x, y, z, 0, ())
plplot.plend()


# vim: expandtab shiftwidth=4:
