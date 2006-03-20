# $Id: ballPlot.py,v 1.6 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting spheres with pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "vtk"
else:
    ren_mod = sys.argv[1]

import random

# set up some data to plot
from Numeric import *

# the three axes in space
# this will give us 10 particles (_not_ 1000)
x = arange(10, typecode=Float)
y = arange(10, typecode=Float)
z = arange(10, typecode=Float)

# 3D position information
posArray = []
for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            posArray.append( (x[i], y[j], z[k]) )

# radius information
random.seed()
radiiArray = zeros(len(x)*len(y)*len(z), typecode=Float)
for i in range(len(x)*len(y)*len(z)):
    radiiArray[i] = random.random()*0.8

# tag information
random.seed()
tagsArray = zeros(len(x)*len(y)*len(z), typecode=Int)
for i in range(len(x)*len(y)*len(z)):
    tagsArray[i] = int(random.random()*10)

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
if ren_mod == "vtk":
    from esys.pyvisi.renderers.vtk import *       # vtk
elif ren_mod == "povray":
    from esys.pyvisi.renderers.povray import *       # povray
else:
    raise ValueError, "Unknown renderer module"

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create a BallPlot object
plot = BallPlot(scene)

# add some helpful info to the plot
plot.title = 'Example ball plot'

# assign some data to the plot
# one way of doing it
# (tags indirectly determine colour of the spheres in the plot)
plot.setData(points=posArray, radii=radiiArray, tags=tagsArray)

# render the scene
scene.render(pause=True, interactive=True)

# without specifying a tags array input
plot.setData(points=posArray, radii=radiiArray)
# render the scene
scene.render(pause=True, interactive=True)

# another way loading an old style-vtk file
plot.setData(fname="cp_test_0.vtk", 
	format="vtk", 
	radii="radius", 
	tags="particleTag")

# render the scene to screen
scene.render(pause=True, interactive=True)

# another way loading a vtk xml file
plot.setData(fname="cp_test_0.xml", 
	format="vtk-xml", 
	radii="radius", 
	tags="particleTag")

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene out to file
scene.save(fname="ballPlot.png", format="png")

# vim: expandtab shiftwidth=4:

