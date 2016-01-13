# Copyright (C) 2004 Paul Cochrane
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# $Id: basicExample.py,v 1.2 2005/01/11 05:43:17 paultcochrane Exp $

## @file basicExample.py

"""
Basic example of pyvisi usage.  

Will hopefully help me write a decent interface.
"""

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

