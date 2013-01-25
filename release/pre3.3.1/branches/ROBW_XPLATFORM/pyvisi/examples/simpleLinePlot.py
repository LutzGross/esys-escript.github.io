# $Id: simpleLinePlot.py,v 1.3 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting lines with pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "vtk"
else:
    ren_mod = sys.argv[1]

# set up some data to plot
from Numeric import *

x = arange(10, typecode=Float)
y = x**2

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
if ren_mod == "gnuplot":
    from esys.pyvisi.renderers.gnuplot import *   # gnuplot
elif ren_mod == "vtk":
    from esys.pyvisi.renderers.vtk import *       # vtk
elif ren_mod == "plplot":
    from esys.pyvisi.renderers.plplot import *    # plplot
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
plot.title = 'Example 2D line plot'
plot.xlabel = 'x'
plot.ylabel = 'x^2'

plot.linestyle = 'lines'

# assign some data to the plot
plot.setData(x, y)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene out to file
scene.save(fname="simpleLinePlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

