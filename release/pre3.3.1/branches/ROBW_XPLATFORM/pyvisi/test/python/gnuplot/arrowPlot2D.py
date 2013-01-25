# $Id: arrowPlot2D.py,v 1.2 2005/05/24 01:46:43 paultcochrane Exp $

"""
Example of plotting a vector field with pyvisi 

This example uses 2D array inputs, which is sometimes easier for users.
"""

# set up some data to plot
from Numeric import *

dim = 10

# initialise the positions of the vectors
x = zeros((dim,dim), typecode=Float)
y = zeros((dim,dim), typecode=Float)

# initialise the vector displacements
# (I may need to rethink how this works in the interface)
dx = zeros((dim,dim), typecode=Float)
dy = zeros((dim,dim), typecode=Float)

# set the positions randomly, and set the displacements to some smaller
# random number but of mean zero instead of distributed between 0 and 1
import random
random.seed()
for i in range(dim):
    for j in range(dim):
        x[i,j] = random.random()
        y[i,j] = random.random()
        dx[i,j] = (random.random()-0.5)/5.0
        dy[i,j] = (random.random()-0.5)/5.0

#### original gnuplot code

import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example 2D arrow/quiver/vector field plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('y')

# set up the data
# first flatten the x,y,dx,dy arrays so can do this
xx = x.flat
yy = y.flat
dxx = dx.flat
dyy = dy.flat
_data = Gnuplot.Data(xx, yy, dxx, dyy, with='vector')

# plot it
_gnuplot.plot(_data)

# set up to save to file
_gnuplot('set terminal png')
_gnuplot('set output \"arrowPlot2D.png\"')

# save it
_gnuplot.plot(_data)

raw_input('Press enter to continue...\n')

