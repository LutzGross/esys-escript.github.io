# Import the necessary modules.
from esys.pyvisi import Scene, DataCollector, Map, Camera, Velocity, Legend 
from esys.pyvisi import Movie, LocalPosition
from esys.pyvisi.constant import *
import os

PYVISI_EXAMPLE_MESHES_PATH = "data_meshes"
PYVISI_EXAMPLE_IMAGES_PATH = "data_sample_images"
X_SIZE = 800
Y_SIZE = 800

SCALAR_FIELD_POINT_DATA = "temp"
FILE_2D = "tempvel-"
IMAGE_NAME = "movie"
JPG_RENDERER = Renderer.ONLINE_JPG

# Create a Scene.
s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, 
        y_size = Y_SIZE)

# Create a DataCollector reading from a XML file. 
dc1 = DataCollector(source = Source.XML)
dc1.setActiveScalar(scalar = SCALAR_FIELD_POINT_DATA)

# Create a Map.
m1 = Map(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, cell_to_point = False,
        outline = True)

# Create a Velocity.
vopc1 = Velocity(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.VECTOR, 
        arrow = Arrow.TWO_D, lut = Lut.COLOR, cell_to_point = False, 
        outline = True)
vopc1.setScaleFactor(scale_factor = 0.07)
vopc1.setRatio(ratio = 8)
vopc1.setColor(color = Color.BLACK)

# Create a scalar Legend.
sb = Legend(scene = s, data_collector= dc1, viewport = Viewport.SOUTH_WEST,
        lut = Lut.COLOR, legend = LegendType.SCALAR)
sb.setOrientationToHorizontal()
sb.setPosition(LocalPosition(85,5))
sb.setScalarRange(0, 1) 
sb.setLabelColor(color = Color.BLACK)
sb.setTitleColor(Color.BLACK)
sb.setTitle("Temperature")

# Create a Camera.
cam1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)

# Create a movie.
mov = Movie()

# Read in one file one after another and render the object. 
for i in range(938, 949):
    dc1.setFileName(file_name =  os.path.join(PYVISI_EXAMPLE_MESHES_PATH, \
	        FILE_2D + "%06d.vtu") % i)

    s.render(image_name = os.path.join(PYVISI_EXAMPLE_IMAGES_PATH, \
            IMAGE_NAME + "%06d.jpg") % i)

# Generate the movie from the rendered images.
mov.makeMovie(input_directory = PYVISI_EXAMPLE_IMAGES_PATH, 
        first_image = IMAGE_NAME + "000938.jpg", 
        last_image = IMAGE_NAME + "000949.jpg", 
        movie = os.path.join(PYVISI_EXAMPLE_IMAGES_PATH, "movie.mpg"))

