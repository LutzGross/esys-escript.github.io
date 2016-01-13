# $Id: seismicOffsetPlotExample2.py,v 1.2 2005/03/22 06:22:36 paultcochrane Exp $ 
"""
Example of plotting multiple curves offset from each other with pyvisi 

This is an example with simulated seismic data, and is a larger dataset
than seismicOffsetPlotExample.py

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


# set up some data to plot
from Numeric import *

# read in the data (being fortunate we know how much data there is)
fp = open('waves1d.dat')
t = zeros((1000), typecode=Float)
x = zeros((102), typecode=Float)
data = zeros((1000,102), typecode=Float)
for i in range(1000):
    for j in range(102):
        line = fp.readline()
        arr = line.split()
        t[i] = float(arr[0])
        x[j] = float(arr[1])
        data[i,j] = float(arr[2])
fp.close()

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
#from esys.pyvisi.utils import *   # pyvisi specific utils
# import the objects to render the scene using the specific renderer
from esys.pyvisi.renderers.gnuplot import *   # gnuplot
#from esys.pyvisi.renderers.vtk import *       # vtk

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create an OffsetPlot object
plot = OffsetPlot(scene)

# add some helpful info to the plot
plot.title = 'Sample seismic data - waves1d.dat'
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
scene.save(fname="seismicOffsetPlotExample2.png", format=PngImage())

# save as postscript
plot.setData(t, data)
                             # have to do this now because we've already
                             # save()d the scene.  This requirement will
                             # be removed in the future.
scene.save(fname="seismicOffsetPlotExample2.ps", format=PsImage())

# vim: expandtab shiftwidth=4:
