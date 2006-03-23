# $Id: scatterPlotExample.py,v 1.4 2005/03/07 04:16:19 paultcochrane Exp $

"""
Example of a scatter plot in pyvisi 

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
import random

x = arange(30, typecode=Float)
y = arange(30, typecode=Float)

# make the data a bit more scatter-like by using random numbers
random.seed()
for i in range(len(x)):
    x[i] = random.random()
    y[i] = random.random()

# plot it using one of the three methods
if method == 'pyvisi':

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
    
    # create a ScatterPlot object
    plot = ScatterPlot(scene)
    
    # add some helpful info to the plot
    plot.title = 'Example 2D scatter plot'
    plot.xlabel = 'x'
    plot.ylabel = 'y'

    # assign some data to the plot
    plot.setData(x, y)

    # render the scene to screen
    scene.render(pause=True, interactive=True)

    # save the scene out to file
    plot.setData(x, y)  # have to do this now because we've already
                        # render()ed the scene.  This requirement will be
                        # removed in the future
    scene.save(fname="scatterPlotExample.png", format=PngImage())
    plot.setData(x, y)  # have to do this now because we've already save()d
                        # the scene.  This requirement will be removed in
                        # the future.
    scene.save(fname="scatterPlotExample.ps", format=PsImage())

elif method == 'gnuplot':
    #### original gnuplot code
    
    import Gnuplot

    # set the plot up
    _gnuplot = Gnuplot.Gnuplot()
    _gnuplot.title('Example 2D scatter plot')
    _gnuplot.xlabel('x')
    _gnuplot.ylabel('y')

    # set up the data
    _data = Gnuplot.Data(x, y, with='points pointtype 2')

    # plot it
    _gnuplot.plot(_data)

    # set up to save to file
    _gnuplot('set terminal png')
    _gnuplot('set output \"scatterPlotExample.png\"')

    # save it
    _gnuplot.plot(_data)

    raw_input('Press enter to continue...\n')

elif method == 'vtk':
    #### original vtk code
    print "vtk scatter plotting not yet implemented"

else:
    print "Eeek!  What plotting method am I supposed to use???"

# vim: expandtab shiftwidth=4:

