from scene import Scene
from datacollector import DataCollector
from map import MapOnPlane
from camera import Camera
from geo import Transform

s = Scene(renderer = "vtk_jpeg", x_size = 800, y_size = 600)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
tf = Transform()
tf.yzPlane(offset = 1)
mop = MapOnPlane(scene = s, data_collector = dc, transform = tf)

# Perform offline rendering  and saving the rendered object as an image.
s.saveImage("save_image.jpg")

