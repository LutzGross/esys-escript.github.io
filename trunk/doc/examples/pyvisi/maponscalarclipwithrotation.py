
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
from esys.escript import getMPIRankWorld
from esys.pyvisi import Scene, DataCollector, MapOnScalarClipWithRotation
from esys.pyvisi import Camera
from esys.pyvisi.constant import *
import os

PYVISI_EXAMPLE_MESHES_PATH = "data_meshes"
PYVISI_EXAMPLE_IMAGES_PATH = "images_out"
if not os.path.isdir(PYVISI_EXAMPLE_IMAGES_PATH) and getMPIRankWorld()==0: os.mkdir(PYVISI_EXAMPLE_IMAGES_PATH)

X_SIZE = 800
Y_SIZE = 800

FILE_2D = "without_st.0700.xml"
IMAGE_NAME = "maponscalarclipwithrotation.jpg"
JPG_RENDERER = Renderer.OFFLINE_JPG # change to Renderer.ONLINE_JPG to interact with visualiztion window

# Create a Scene.
s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
        y_size = Y_SIZE)

# Create a DataCollector reading from a XML file.
dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = os.path.join(PYVISI_EXAMPLE_MESHES_PATH, FILE_2D))

# Create a MapOnScalarClipWithRotation.
m1 = MapOnScalarClipWithRotation(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, cell_to_point = False)
m1.setAngle(220)
m1.setResolution(50)

# Create a Camera.
c2 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
c2.isometricView()
c2.elevation(-30)

# Render the object.
s.render(image_name = os.path.join(PYVISI_EXAMPLE_IMAGES_PATH, IMAGE_NAME))
