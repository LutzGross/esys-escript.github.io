# $Id: scatterPlot.py,v 1.1 2005/11/08 07:34:47 paultcochrane Exp $

"""
Example of a scatter plot
"""

# set up some data to plot
from Numeric import *
import random

x = arange(30, typecode=Float)
y = arange(30, typecode=Float)

# make the data a bit more scatter-like by using random numbers
random.seed()
for i in range(len(x)):
    x[i] = random.random()
    y[i] = random.random()

#### original gnuplot code

import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example 2D scatter plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('y')

# set up the data
_data = Gnuplot.Data(x, y, with='points pointtype 2')

# plot it
_gnuplot.plot(_data)

# set up to save to file
_gnuplot('set terminal png')
_gnuplot('set output \"scatterPlot.png\"')

# save it
_gnuplot.plot(_data)

raw_input('Press enter to continue...\n')

# vim: expandtab shiftwidth=4:

