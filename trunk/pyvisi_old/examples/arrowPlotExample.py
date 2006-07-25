# $Id: arrowPlotExample.py,v 1.3 2005/05/24 01:28:19 paultcochrane Exp $

"""
Example of plotting a vector field with pyvisi 

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


# what plotting method are we using?
method = 'pyvisi'

# set up some data to plot
from Numeric import *

# the positions of the vectors
x = arange(10, typecode=Float)
y = arange(10, typecode=Float)

# the vector displacements
# (I may need to rethink how this works in the interface)
dx = arange(10, typecode=Float)
dy = arange(10, typecode=Float)

# set the positions randomly, and set the displacements to be the square of
# the positions
import random
random.seed()
for i in range(len(x)):
    x[i] = random.random()
    y[i] = random.random()
    dx[i] = x[i]*x[i]
    dy[i] = y[i]*y[i]

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
    scene.save(fname="arrowPlotExample.png", format=PngImage())

elif method == 'vtk':
    #### original vtk code

    print "Sorry, the vtk interface hasn't been written yet."
elif method == 'plplot':
    #### original plplot code

    print "Sorry, the plplot interface hasn't been written yet."
else:
    print "Eeek!  What plotting method am I supposed to use???"

# vim: expandtab shiftwidth=4:

