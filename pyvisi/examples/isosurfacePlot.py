# $Id: isosurfacePlot.py,v 1.4 2005/11/08 08:23:45 paultcochrane Exp $

"""
Example of plotting a set of isosurfaces with pyvisi 

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

# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
# import the objects to render the scene using the specific renderer
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

# create a IsosurfacePlot object
plot = IsosurfacePlot(scene)

# add some helpful info to the plot
plot.title = 'Example isosurface plot'
plot.xlabel = 'x'
plot.ylabel = 'y'
plot.zlabel = 'z'

# plot data defined in a vtk file
plot.setData(fname='temp-0500.vtk', format='vtk-xml')

scene.render(pause=True, interactive=True)

# save the plot
scene.save(fname="isosurfacePlot.png", format="png")

# vim: expandtab shiftwidth=4:

