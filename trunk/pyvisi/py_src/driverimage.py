from scene import Scene
from image import Image

s = Scene(renderer = "vtk_online", x_size = 800, y_size = 600)
# Display an image.
i = Image(scene = s, format = "jpeg")
i.setFileName("../test/python/data_data/Flinders_eval.jpg" )
s.render()
