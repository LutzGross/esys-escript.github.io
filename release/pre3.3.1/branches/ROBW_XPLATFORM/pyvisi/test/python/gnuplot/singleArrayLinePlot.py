# $Id: singleArrayLinePlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

"""
Example of plotting a curve using only one input array with pyvisi 
"""

# set up some data to plot
from Numeric import *

x = arange(0,2*pi,0.1, typecode=Float)
y = sin(x)

#### original gnuplot code
    
import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example 2D plot')
_gnuplot.xlabel('index')
_gnuplot.ylabel('y')

# set up the data
_data = Gnuplot.Data(y, with='lines')

# plot it
_gnuplot.plot(_data)

# save it to file
_gnuplot('set terminal png')
_gnuplot('set output "singleArrayLinePlot.png"')
_gnuplot.plot(_data)

raw_input('Press enter to continue...\n')


