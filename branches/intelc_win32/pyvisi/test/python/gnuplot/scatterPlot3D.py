# $Id: scatterPlot3D.py,v 1.1 2005/11/08 07:39:37 paultcochrane Exp $

"""
Example of plotting 3D scatter plots 

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
import random

# the x and y axes
x = arange(-2, 2, 0.2, typecode=Float)
y = arange(-3, 3, 0.2, typecode=Float)

# make the data a bit more scatter-like by using random numbers
random.seed()
for i in range(len(x)):
    x[i] = x[i]*random.random()

for i in range(len(y)):
    y[i] = y[i]*random.random()

# pick some interesting function to generate the data in the third dimension
# this is the one used in the matlab docs: z = x*exp(-x^2-y^2)
z = zeros((len(x),len(y)), typecode=Float)

# boy do *I* feel old fashioned writing it this way
# surely there's another way to do it: - something to do later
for i in range(len(x)):
    for j in range(len(y)):
	z[i,j] = exp(-x[i]*x[i] - y[j]*y[j])

#### original gnuplot code

import Gnuplot

# set the plot up
_gnuplot = Gnuplot.Gnuplot()
_gnuplot.title('Example 3D scatter plot')
_gnuplot.xlabel('x')
_gnuplot.ylabel('y')
_gnuplot("set zlabel \'z\'")

# scatter plot, so set the data style to points
_gnuplot('set data style points')

# set up the data
_data = Gnuplot.GridData(z,x,y, binary=1)

_gnuplot.splot(_data)

raw_input('Press enter to continue...')

# vim: expandtab shiftwidth=4: