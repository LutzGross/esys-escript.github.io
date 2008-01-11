"""
Author: John Ngui, john.ngui@uq.edu.au
"""

# Import the necessary modules.
from esys.pyvisi import Scene, DataCollector, Map, Camera, Legend, Contour
from esys.pyvisi import LocalPosition
from esys.pyvisi.constant import *
import os

PYVISI_EXAMPLE_MESHES_PATH = "data_meshes"
PYVISI_EXAMPLE_IMAGES_PATH = "data_sample_images"
X_SIZE = 500
Y_SIZE = 500

FILE_NAME = "phi_talus_lava.0099.vtu"
IMAGE_NAME = "legend.jpg"
JPG_RENDERER = Renderer.ONLINE_JPG

# Create a Scene.
s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
        y_size = Y_SIZE)

# Create a DataCollector reading from a XML file. 
dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = os.path.join(PYVISI_EXAMPLE_MESHES_PATH, \
        FILE_NAME))

# Create a Contour.
ctr1 = Contour(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, cell_to_point = False,
        outline = True)
ctr1.generateContours(8, -4.963259999999998, 86.277230000000003)
ctr1.setScalarRange(-4.963259999999998, 86.277230000000003)

# Create a scalar Legend.
lg1 =Legend(scene = s, data_collector= dc1, viewport = Viewport.SOUTH_WEST,
        lut = Lut.COLOR, legend = LegendType.SCALAR)
lg1.setOrientationToHorizontal()
lg1.setScalarRange(-4.963259999999998, 86.277230000000003)
lg1.setTitle(title = "Scalar Bar")
lg1.setPosition(LocalPosition(50, 5))

# Create a Camera.
cam1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)

# Render the object. 
s.render(image_name = os.path.join(PYVISI_EXAMPLE_IMAGES_PATH, 
        IMAGE_NAME))
