# $Id: surfacePlot.py,v 1.5 2005/11/10 08:38:06 paultcochrane Exp $

"""
Example of plotting surfaces with pyvisi 

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


import sys
numArgs = len(sys.argv)
if numArgs == 1:
    ren_mod = "vtk"
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
from esys.pyvisi import *
# import the gnuplot overrides of the interface
if ren_mod == "gnuplot":
    from esys.pyvisi.renderers.gnuplot import *
elif ren_mod == "vtk":
    from esys.pyvisi.renderers.vtk import *
elif ren_mod == "plplot":
    from esys.pyvisi.renderers.plplot import *
else:
    raise ValueError, "Unknown renderer module"

# define a scene object
# a Scene is a container for all of the kinds of things you want to put
# into your plot, for instance, images, meshes, arrow/vector/quiver
# plots, contour plots, spheres etc.
scene = Scene()

# create a SurfacePlot object
plot = SurfacePlot(scene)

# add some helpful info to the plot
plot.title = 'Example surface plot'
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
scene.save(fname="surfacePlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:
