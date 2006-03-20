# $Id: arrowPlotExample2.py,v 1.2 2005/05/24 01:30:36 paultcochrane Exp $

"""
Example of plotting a vector field with pyvisi 

This example uses 2D array inputs, which is sometimes easier for users.
"""

# what plotting method are we using?
method = 'pyvisi'

# set up some data to plot
from Numeric import *

dim = 10

# initialise the positions of the vectors
x = zeros((dim,dim), typecode=Float)
y = zeros((dim,dim), typecode=Float)

# initialise the vector displacements
# (I may need to rethink how this works in the interface)
dx = zeros((dim,dim), typecode=Float)
dy = zeros((dim,dim), typecode=Float)

# set the positions randomly, and set the displacements to some smaller
# random number but of mean zero instead of distributed between 0 and 1
import random
random.seed()
for i in range(dim):
    for j in range(dim):
        x[i,j] = random.random()
        y[i,j] = random.random()
        dx[i,j] = (random.random()-0.5)/5.0
        dy[i,j] = (random.random()-0.5)/5.0

# plot it using one of the three methods
if method == 'pyvisi':

    # example code for how a user would write a script in pyvisi
    from esys.pyvisi import *          # base level visualisation stuff
    #from esys.pyvisi.utils import *   # pyvisi specific utils
    # import the objects to render the scene using the specific renderer
    from esys.pyvisi.renderers.gnuplot import *   # gnuplot
    
    # define the scene object
    # a Scene is a container for all of the kinds of things you want to put 
    # into your plot for instance, images, meshes, arrow/vector/quiver plots, 
    # contour plots, spheres etc.
    scene = Scene()
    
    # create a LinePlot object
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
    plot.setData(x, y, dx, dy) # have to do this now because we've already
                               # render()ed the scene.  This requirement
                               # will be removed in the future
    scene.save(fname="arrowPlotExample2.png", format=PngImage())

elif method == 'vtk':
    #### original vtk code

    print "Sorry, the vtk interface hasn't been written yet."
elif method == 'plplot':
    #### original plplot code

    print "Sorry, the plplot interface hasn't been written yet."
else:
    print "Eeek!  What plotting method am I supposed to use???"

# vim: expandtab shiftwidth=4:

