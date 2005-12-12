# $Id: scatterPlot.py,v 1.3 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of a scatter plot in pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "gnuplot"
else:
    ren_mod = sys.argv[1]

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

# example code for how a user would write a script in pyvisi
from pyvisi import *          # base level visualisation stuff
#from pyvisi.utils import *   # pyvisi specific utils
# import the objects to render the scene using the specific renderer
if ren_mod == "gnuplot":
    from pyvisi.renderers.gnuplot import *   # gnuplot
elif ren_mod == "vtk":
    from pyvisi.renderers.vtk import *       # vtk
else:
    raise ValueError, "Unknown renderer module"

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create a ScatterPlot object
plot = ScatterPlot(scene)

# add some helpful info to the plot
plot.title = 'Example 2D scatter plot'
plot.xlabel = 'x'
plot.ylabel = 'y'

# assign some data to the plot
plot.setData(x, y)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene out to file
scene.save(fname="scatterPlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

