# $Id: singleArrayLinePlot.py,v 1.3 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting a curve using only one input array with pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "vtk"
else:
    ren_mod = sys.argv[1]

# set up some data to plot
from Numeric import *

x = arange(0,2*pi,0.1, typecode=Float)
y = sin(x)

# example code for how a user would write a script in pyvisi
from pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
if ren_mod == "gnuplot":
    from pyvisi.renderers.gnuplot import *   # gnuplot
elif ren_mod == "vtk":
    from pyvisi.renderers.vtk import *       # vtk
elif ren_mod == "plplot":
    from pyvisi.renderers.plplot import *    # plplot
else:
    raise ValueError, "Unknown renderer module"

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create a LinePlot object
plot = LinePlot(scene)

# add some helpful info to the plot
plot.title = 'Example 2D plot'
plot.xlabel = 'index'
plot.ylabel = 'y'

plot.linestyle = 'lines'

# assign some data to the plot
plot.setData(y)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene to file
scene.save(fname="singleArrayLinePlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

