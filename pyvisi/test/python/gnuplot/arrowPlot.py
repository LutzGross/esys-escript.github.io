# $Id: arrowPlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

"""
Example of plotting a vector field with pyvisi 


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

# the positions of the vectors
x = arange(10, typecode=Float)
y = arange(10, typecode=Float)

# the vector displacements
# (I may need to rethink how this works in the interface)
dx = arange(10, typecode=Float)
dy = arange(10, typecode=Float)

# set the positions randomly, and set the displacements to be the square of
# the positions
import random
random.seed()
for i in range(len(x)):
    x[i] = random.random()
    y[i] = random.random()
    dx[i] = x[i]*x[i]
    dy[i] = y[i]*y[i]

#### original gnuplot code
    
import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example 2D arrow/quiver/vector field plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('y')

# set up the data
_data = Gnuplot.Data(x, y, dx, dy, with='vector')

# plot it
_gnuplot.plot(_data)

# set up to save to file
_gnuplot('set terminal png')
_gnuplot('set output \"arrowPlotExample.png\"')

# save it
_gnuplot.plot(_data)

raw_input('Press enter to continue...\n')


