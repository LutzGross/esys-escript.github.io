# $Id: offsetLinePlot.py,v 1.4 2005/11/08 08:23:45 paultcochrane Exp $ 
"""
Example of plotting multiple curves offset from each other with pyvisi 

This is especially handy for people plotting seismic data

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

x = arange(0,2*pi,0.01, typecode=Float)
y1 = sin(x)
y2 = cos(x)
y3 = cos(x)**2
y4 = sin(2*x)
y5 = cos(3*x)
y6 = sin(20*x)

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
#from esys.pyvisi.utils import *   # pyvisi specific utils
# import the objects to render the scene using the specific renderer
if ren_mod == "gnuplot":
    from esys.pyvisi.renderers.gnuplot import *   # gnuplot
elif ren_mod == "vtk":
    from esys.pyvisi.renderers.vtk import *       # vtk
elif ren_mod == "plplot":
    from esys.pyvisi.renderers.plplot import *    # plplot
else:
    raise ValueError, "Unknown renderer module"

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# create a LinePlot object
plot = LinePlot(scene)

# add some helpful info to the plot
plot.title = 'Example 2D plot with offsets'
plot.xlabel = 'x'
plot.ylabel = 'y'

plot.linestyle = 'lines'

# assign some data to the plot
plot.setData(x, y1, y2, y3, y4, y5, y6, offset=True)

# render the scene to screen
scene.render(pause=True, interactive=True)

# save the scene to file
scene.save(fname="offsetLinePlot.png", format=PngImage())

# vim: expandtab shiftwidth=4:

