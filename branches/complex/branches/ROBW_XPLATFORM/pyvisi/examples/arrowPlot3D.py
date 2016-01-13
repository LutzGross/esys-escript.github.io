# $Id: arrowPlot3D.py,v 1.6 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting a 3D vector field with pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "vtk"
else:
    ren_mod = sys.argv[1]

# set up some data to plot
from Numeric import *

dim = 10

# initialise the positions of the vectors
x = zeros((dim,dim), typecode=Float)
y = zeros((dim,dim), typecode=Float)
z = zeros((dim,dim), typecode=Float)

# initialise the vector displacements
# (I may need to rethink how this works in the interface)
dx = zeros((dim,dim), typecode=Float)
dy = zeros((dim,dim), typecode=Float)
dz = zeros((dim,dim), typecode=Float)

# set the positions randomly, and set the displacements to some smaller
# random number but of mean zero instead of distributed between 0 and 1
import random
random.seed()
for i in range(dim):
    for j in range(dim):
        x[i,j] = random.random()
        y[i,j] = random.random()
        z[i,j] = random.random()
        dx[i,j] = (random.random()-0.5)/5.0
        dy[i,j] = (random.random()-0.5)/5.0
        dz[i,j] = (random.random()-0.5)/5.0

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
#from esys.pyvisi.renderers.gnuplot import *   # gnuplot
if ren_mod == "vtk":
    from esys.pyvisi.renderers.vtk import *       # vtk
elif ren_mod == "povray":
    from esys.pyvisi.renderers.povray import *    # povray
else:
    raise ValueError, "Unknown renderer module"

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create a ArrowPlot3D object
plot = ArrowPlot3D(scene)

# add some helpful info to the plot
plot.title = 'Example 3D arrow/quiver/vector field plot'
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.zlabel = 'z'

# assign some data to the plot
plot.setData(x, y, z, dx, dy, dz)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene out to file
scene.save(fname="arrowPlot3D.png", format=PngImage())

# plot data defined in a vtk file
plot.setData(fname='vel-0500.vtk', format='vtk-xml')

scene.render(pause=True, interactive=True)

# save this plot too
scene.save(fname="arrowPlot3Dvtk.png", format="png")

# vim: expandtab shiftwidth=4:

