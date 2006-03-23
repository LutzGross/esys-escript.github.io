# $Id: singleArrayLinePlot.py,v 1.1 2005/05/05 01:56:51 paultcochrane Exp $

"""
Example of plotting a curve using only one input array with pyvisi 

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


