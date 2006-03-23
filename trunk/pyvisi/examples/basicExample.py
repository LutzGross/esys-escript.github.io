"""
Basic example of pyvisi usage.  

Will hopefully help me write a decent interface.

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


# example code for how a user would write a script in pyvisi
from esys.pyvisi import *          # base level visualisation stuff
#from esys.pyvisi.utils import *   # pyvisi specific utils
# import the objects to render the scene using vtk
from esys.pyvisi.renderers.vtk import * 

# these things are just here to make data to plot, not all of which are used
from ESyS import *
import Finley

# now make some data of some kind
mesh = Finley.Brick(3,5,7)  # a Finley mesh
vectorData = mesh.Nodes().getX()  # get vector data from the mesh nodes

# define the scene object
# a Scene is a container for all of the kinds of things you want to put 
# into your plot for instance, images, meshes, arrow/vector/quiver plots, 
# contour plots, spheres etc.
scene = Scene()

# define a camera object.  There will need to be one camera per scene.
camera = Camera()

# add the camera to the scene
scene.add(camera)

# create an ArrowPlot object
#plot = ArrowPlot()

# add the plot to the scene
#scene.add(plot)

# assign some data to the plot
#plot.setData(vectorData)

# create an Image object
img = Image(file="ranges.jpg",format="jpeg")

# add the image to the scene
scene.add(image)

# render the scene, outputing the data to a jpeg file
scene.render(file="example.jpg",format="jpeg")
# if we are just using vtk, this should head out to a file, or an opengl window
# if we are using mayavi with vtk, then this will be in a mayavi window

# saving a scene could also be written as (handy for PBS jobs?)
# could try working out the format from the filename extension
scene.save(file="example.jpg", format="jpeg") 

