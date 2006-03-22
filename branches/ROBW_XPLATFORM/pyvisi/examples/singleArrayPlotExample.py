# $Id: singleArrayPlotExample.py,v 1.12 2005/05/24 01:44:15 paultcochrane Exp $

"""
Example of plotting a curve using only one input array with pyvisi 
"""

# set up some data to plot
from Numeric import *

x = arange(0,2*pi,0.1, typecode=Float)
y = sin(x)

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
plot.xlabel = 'index'
plot.ylabel = 'y'

plot.linestyle = 'lines'

# assign some data to the plot
plot.setData(y)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene to file
plot.setData(y)  # have to do this now because we've already render()ed
                 # the scene.  This requirement will be removed in the
                 # future
scene.save(fname="singleArrayLinePlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

