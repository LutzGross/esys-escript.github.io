from scene import Scene
from datacollector import DataCollector
from camera import Camera
from contour import Contour

s = Scene(renderer = "vtk_online", x_size = 800, y_size = 600)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -30)

# Creates contours.
c = Contour(scene = s, data_collector = dc)
# Specify the number of contours and the range.
c.generateValues(number_contours = 5, min_range = 0.0, max_range = 1.2)
s.render()

