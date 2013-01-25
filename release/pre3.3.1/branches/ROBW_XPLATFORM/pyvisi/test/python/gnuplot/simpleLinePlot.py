# $Id: simpleLinePlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

"""
Example of plotting lines with pyvisi

This is the original code used to develop the gnuplot renderer module
"""

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

