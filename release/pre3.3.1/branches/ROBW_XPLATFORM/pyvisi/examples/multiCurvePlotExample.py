# $Id: multiCurvePlotExample.py,v 1.19 2005/05/24 01:36:25 paultcochrane Exp $

"""
Example of plotting multiple curves with pyvisi 
"""

# set up some data to plot
from Numeric import *

x = arange(0, 2*pi, 0.1, typecode=Float)
y1 = sin(x)
y2 = cos(x)
y3 = cos(x)**2

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
from esys.pyvisi.renderers.gnuplot import *   # gnuplot
#from esys.pyvisi.renderers.vtk import *       # vtk
#from esys.pyvisi.renderers.plplot import *    # plplot

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create a LinePlot object
plot = LinePlot(scene)

# add some helpful info to the plot
plot.title = 'Example 2D plot'
plot.xlabel = 'x'
plot.ylabel = 'y'

plot.linestyle = 'lines'

# assign some data to the plot
plot.setData(x, y1, y2, y3)

# render the scene to screen
scene.render(pause=True,interactive=True)

# save the scene to file
plot.setData(x, y1, y2, y3)  # have to do this now because we've already
                             # render()ed the scene.  This requirement
                             # will be removed in the future.
scene.save(fname="multiCurveLinePlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

