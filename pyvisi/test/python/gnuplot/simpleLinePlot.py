# $Id: simpleLinePlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

"""
Example of plotting lines with pyvisi

This is the original code used to develop the gnuplot renderer module

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

x = arange(10, typecode=Float)
y = x**2

import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example 2D plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('x^2')

# set up the data
_data = Gnuplot.Data(x, y, with='lines')

# plot it
_gnuplot.plot(_data)

# set up to save to file
_gnuplot('set terminal png')
_gnuplot('set output \"simpleLinePlot.png\"')

# save it
_gnuplot.plot(_data)

raw_input('Press enter to continue...\n')

