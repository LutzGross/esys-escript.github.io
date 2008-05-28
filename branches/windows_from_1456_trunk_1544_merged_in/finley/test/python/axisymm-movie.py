from esys.pyvisi import Scene, DataCollector, Map, Camera, Velocity, Legend 
from esys.pyvisi import Movie, LocalPosition
from esys.pyvisi.constant import *
import os

X_SIZE = 800
Y_SIZE = 600


JPG_RENDERER = Renderer.ONLINE_JPG
s = Scene(renderer = JPG_RENDERER, num_viewport = 1, x_size = X_SIZE, y_size = Y_SIZE)

# Create a DataCollector reading from a XML file. 
dc1 = DataCollector(source = Source.XML)
dc1.setActiveScalar(scalar = "p")
dc1.setActiveVector(vector = "U")

# Create a Map.
m1 = Map(scene = s, data_collector = dc1, lut = Lut.COLOR, cell_to_point = False, outline = True)

vopc1 = Velocity(scene = s, data_collector = dc1,
        color_mode = ColorMode.VECTOR,
        arrow = Arrow.THREE_D, lut = Lut.COLOR, cell_to_point = False,
        outline = False)
vopc1.setScaleFactor(scale_factor = 0.04)


# Create a Camera.
cam1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)

# Create a movie.
mov = Movie()
#lst = []

# Read in one file one after another and render the object. 
images=[]
for i in range(0,50):
    dc1.setFileName("u.%d.xml"%i)
    image="frame.%06d.jpg"%i
    s.render(image_name = os.path.join(".",image))
    images.append(image)

mov.imageList(input_directory = ".", image_list = images)

# Generate the movie.
mov.makeMovie(os.path.join(".", "movie.mpg"))

