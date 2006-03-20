# $Id: multiCurveLinePlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

"""
Example of plotting multiple curves with pyvisi 
"""

# set up some data to plot
from Numeric import *

x = arange(0, 2*pi, 0.1, typecode=Float)
y1 = sin(x)
y2 = cos(x)
y3 = cos(x)**2

#### original gnuplot code

import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example 2D plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('y')

# set up the data
_data1 = Gnuplot.Data(x, y1, with='lines')
_data2 = Gnuplot.Data(x, y2, with='lines')
_data3 = Gnuplot.Data(x, y3, with='lines')

# plot it
_gnuplot.plot(_data1, _data2, _data3)

# save it to file
_gnuplot('set terminal png')
_gnuplot('set output "multiCurveLinePlot.png"')
_gnuplot.plot(_data1, _data2, _data3)

raw_input('Press enter to continue...\n')

