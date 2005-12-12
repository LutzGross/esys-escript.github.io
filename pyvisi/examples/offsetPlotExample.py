# $Id: offsetPlotExample.py,v 1.4 2005/05/24 01:37:50 paultcochrane Exp $ 
"""
Example of plotting multiple curves offset from each other with pyvisi 

This is an example with simulated seismic data, and is a larger dataset
than seismicOffsetPlotExample.py
"""

# set up some data to plot
from Numeric import *

# read in the data (being fortunate we know how much data there is)
fp = open('waves.dat')
t = zeros((100), typecode=Float)
x = zeros((13), typecode=Float)
data = zeros((100,13), typecode=Float)
for i in range(100):
    for j in range(13):
        line = fp.readline()
        arr = line.split()
        t[i] = float(arr[0])
        x[j] = float(arr[1])
        data[i,j] = float(arr[2])
fp.close()

# example code for how a user would write a script in pyvisi
from pyvisi import *          # base level visualisation stuff
#from pyvisi.utils import *   # pyvisi specific utils
# import the objects to render the scene using the specific renderer
#from pyvisi.renderers.gnuplot import *   # gnuplot
from pyvisi.renderers.vtk import *       # vtk

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create an OffsetPlot object
plot = OffsetPlot(scene)

# add some helpful info to the plot
plot.title = 'OffsetPlot example - waves.dat'
plot.xlabel = 't'
plot.ylabel = 'y'

# assign some data to the plot
plot.setData(t, data)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene to file
# save as png
plot.setData(t, data)
                             # have to do this now because we've already
                             # render()ed the scene.  This requirement
                             # will be removed in the future.
scene.save(fname="offsetPlotExample.png", format=PngImage())

# vim: expandtab shiftwidth=4:

