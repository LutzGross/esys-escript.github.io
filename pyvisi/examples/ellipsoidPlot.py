# $Id: ellipsoidPlot.py,v 1.5 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting ellipsoids (useful for visualising tensors) with pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "vtk"
else:
    ren_mod = sys.argv[1]

# set up some data to plot
from Numeric import *

# example code for how a user would write a script in pyvisi
from pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer

if ren_mod == "vtk":
    from pyvisi.renderers.vtk import *       # vtk
elif ren_mod == "povray":
    from pyvisi.renderers.povray import *    # povray
else:
    raise ValueError, "Unknown renderer module"

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create a EllipsoidPlot object
plot = EllipsoidPlot(scene)

# add some helpful info to the plot
plot.title = 'Example ellipsoid plot'

# plot data defined in a vtk file
plot.setData(fname='stress22.vtk', format='vtk-xml')

scene.render(pause=True, interactive=True)

# save the plot
scene.save(fname="ellipsoidPlot.png", format="png")

# vim: expandtab shiftwidth=4:

