# $Id: meshPlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

"""
Example of plotting meshed surfaces with pyvisi 
"""

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
_gnuplot.title('Example mesh plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('y')
_gnuplot('set zlabel \'z\'')

# this is a mesh plot, so...
_gnuplot('set surface')
_gnuplot('set data style lines')

# set up the data
_data = Gnuplot.GridData(z,x,y, binary=1)

_gnuplot.splot(_data)

# set up to save to file
_gnuplot('set terminal png')
_gnuplot('set output \"meshPlot.png\"')

# save it
_gnuplot.splot(_data)

raw_input('Press enter to continue...')

