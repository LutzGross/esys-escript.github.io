from scene import Scene
from datacollector import DataCollector
from camera import Camera
from tensor import Tensor

s = Scene(renderer = "vtk_online", x_size = 800, y_size = 600)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -50)

# Create ellipsoids.
t = Tensor(scene = s, data_collector = dc)
# Change the tensor's default scale factor.
t.setScaleFactor(scale_factor = 0.1)
s.render()

