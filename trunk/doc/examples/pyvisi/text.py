
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: John Ngui, john.ngui@uq.edu.au
"""

# Import the necessary modules.
from esys.pyvisi import Scene, Text2D, LocalPosition
from esys.pyvisi.constant import *
from esys.escript import getMPIRankWorld
import os

PYVISI_EXAMPLE_MESHES_PATH = "data_meshes"
PYVISI_EXAMPLE_IMAGES_PATH = "images_out"
if not os.path.isdir(PYVISI_EXAMPLE_IMAGES_PATH) and getMPIRankWorld()==0: os.mkdir(PYVISI_EXAMPLE_IMAGES_PATH)

X_SIZE = 600
Y_SIZE = 600

IMAGE_NAME = "text.jpg"
JPG_RENDERER = Renderer.OFFLINE_JPG # change to Renderer.ONLINE_JPG to interact with visualiztion window

# Create a Scene.
s = Scene(renderer = JPG_RENDERER, num_viewport = 4, x_size = X_SIZE, 
        y_size = Y_SIZE)

# Create a 2D text for the first viewport.
t1 = Text2D(scene = s, text = "VTK ...", viewport = Viewport.SOUTH_WEST)
t1.setPosition(LocalPosition(x_coor = 20, y_coor = 30))
t1.setColor(color = Color.BLUE)
t1.setFontSize(size = 20)
t1.boldOn()

# Create a 2D text for the third viewport.
t2 = Text2D(scene = s, text = "PYTHON ...\n MESH", 
        viewport = Viewport.NORTH_EAST)
t2.setPosition(LocalPosition(x_coor = 20, y_coor = 30))
t2.setColor(color = Color.PURPLE)
t2.setFontSize(size = 50)
t2.setFontToArial()
t2.shadowOn()

# Render the object.
s.render(image_name = os.path.join(PYVISI_EXAMPLE_IMAGES_PATH, IMAGE_NAME))

