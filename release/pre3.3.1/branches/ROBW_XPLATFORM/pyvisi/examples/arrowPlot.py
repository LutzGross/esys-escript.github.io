# $Id: arrowPlot.py,v 1.4 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting a vector field with pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "vtk"
else:
    ren_mod = sys.argv[1]

# set up some data to plot
from Numeric import *

# the positions of the vectors
x = arange(20, typecode=Float)
y = arange(20, typecode=Float)

# the vector displacements
# (I may need to rethink how this works in the interface)
dx = arange(20, typecode=Float)
dy = arange(20, typecode=Float)

# set the positions randomly, and set the displacements to be the square of
# the positions
import random
random.seed()
for i in range(len(x)):
    x[i] = random.random()
    y[i] = random.random()
    dx[i] = x[i]*x[i]
    dy[i] = y[i]*y[i]

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
if ren_mod == "gnuplot":
    from esys.pyvisi.renderers.gnuplot import *   # gnuplot
elif ren_mod == "vtk":
    from esys.pyvisi.renderers.vtk import *   # vtk
else:
    raise ValueError, "Unknown renderer module"

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create an ArrowPlot object
plot = ArrowPlot(scene)

# add some helpful info to the plot
plot.title = 'Example 2D arrow/quiver/vector field plot'
plot.xlabel = 'x'
plot.ylabel = 'y'

# assign some data to the plot
plot.setData(x, y, dx, dy)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene out to file
scene.save(fname="arrowPlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

