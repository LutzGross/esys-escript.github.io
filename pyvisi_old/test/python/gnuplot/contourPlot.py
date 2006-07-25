# $Id: contourPlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

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

#### original gnuplot code

import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example contour plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('y')

# this is a contour plot, so...
_gnuplot('set contour base')
_gnuplot('set view 0, 0, 1, 1')
_gnuplot('set nosurface')

# set up the data
_data = Gnuplot.GridData(z,x,y, binary=1)

_gnuplot.splot(_data)

# set up to save to file
_gnuplot('set terminal png')
_gnuplot('set output \"contourPlot.png\"')

# save it
_gnuplot.splot(_data)

raw_input('Press enter to continue...')

