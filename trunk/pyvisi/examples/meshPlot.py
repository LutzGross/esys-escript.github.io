# $Id: meshPlot.py,v 1.3 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting meshed surfaces with pyvisi 
"""

import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "gnuplot"
else:
    ren_mod = sys.argv[1]

# set up some data to plot
from Numeric import *

# the x and y axes
x = arange(-2,2,0.2, typecode=Float)
y = arange(-2,3,0.2, typecode=Float)

# pick some interesting function to generate the data in the third dimension
# this is the one used in the matlab docs: z = x*exp(-x^2-y^2)
z = zeros((len(x),len(y)), typecode=Float)

# boy do *I* feel old fashioned writing it this way
# surely there's another way to do it: - something to do later
for i in range(len(x)):
    for j in range(len(y)):
	z[i,j] = x[i]*exp(-x[i]*x[i] - y[j]*y[j])

# import the general pyvisi stuff
from pyvisi import *
if ren_mod == "gnuplot":
    from pyvisi.renderers.gnuplot import *
elif ren_mod == "vtk":
    from pyvisi.renderers.vtk import *
else:
    raise ValueError, "Unknown renderer module"

# define a scene object
# a Scene is a container for all of the kinds of things you want to put
# into your plot, for instance, images, meshes, arrow/vector/quiver
# plots, contour plots, spheres etc.
scene = Scene()

# create a MeshPlot object
plot = MeshPlot(scene)

# add some helpful info to the plot
plot.title = 'Example mesh plot'
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.zlabel = 'z'

# assign the data to the plot
# this version assumes that we have x, then y, then z and that z is 2D
# and that x and y are 1D arrays
plot.setData(x,y,z)
# alternative syntax
#plot.setData(xData=x, yData=y, zData=z)
# or (but more confusing depending upon one's naming conventions)
#plot.setData(x=x, y=y, z=z)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene to file
scene.save(fname="meshPlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:
