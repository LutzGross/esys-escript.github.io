# $Id: simplePlotExample.py,v 1.9 2005/05/24 01:42:44 paultcochrane Exp $

"""
Example of plotting lines with pyvisi 
"""

# set up some data to plot
from Numeric import *

x = arange(10, typecode=Float)
y = x**2

# example code for how a user would write a script in pyvisi
from pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
#from pyvisi.renderers.gnuplot import *   # gnuplot
from pyvisi.renderers.vtk import *       # vtk
#from pyvisi.renderers.plplot import *    # plplot

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
## png
plot.setData(x, y)  # have to do this now because we've already
                    # render()ed the scene, will be removed in the
                    # future
scene.save(fname="simpleLinePlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

