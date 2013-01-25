# $Id: simpleLinePlot.py,v 1.1 2005/05/05 01:57:24 paultcochrane Exp $

"""
Example of plotting lines with pyvisi

This is the original code used to develop the plplot renderer module
"""

# set up some data to plot
from Numeric import *

x = arange(10, typecode=Float)
y = x**2

import plplot

plplot.plsdev("xwin")
plplot.plinit()
plplot.plenv(min(x), max(x), min(y), max(y), 0, 1)
plplot.pllab("x", "x**2", "Example 2D plot")
plplot.plline(x, y)
plplot.plend()

# to save as well, have to set everything up again, and replot
# save as png
plplot.plsdev("png")
plplot.plsfnam("simplePlotExample.png")
plplot.plinit()
plplot.plenv(min(x), max(x), min(y), max(y), 0, 1)
plplot.pllab("x", "x**2", "Example 2D plot")
plplot.plline(x, y)
plplot.plend()

# save as postscript
plplot.plsdev("psc")
plplot.plsfnam("simplePlotExample.ps")
plplot.plinit()
plplot.plenv(min(x), max(x), min(y), max(y), 0, 1)
plplot.pllab("x", "x**2", "Example 2D plot")
plplot.plline(x, y)
plplot.plend()
