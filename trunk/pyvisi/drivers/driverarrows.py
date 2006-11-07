from scene import Scene
from datacollector import DataCollector
from camera import Camera
from arrows import Arrows 

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
#dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")
dc.setFileName(file_name = "../../../data/laurent/slab/slab.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -40)

# Create arrows.
a = Arrows(scene = s, data_collector = dc)
# Change the arrows' default scale_factor
#a.setScaleFactor(scale_factor = 0.25)
a.setScaleFactor(scale_factor = 4000)
s.render()

