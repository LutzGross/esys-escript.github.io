from scene import Scene
from datacollector import DataCollector
from camera import Camera
from isosurface import IsoSurface

s = Scene(renderer = "vtk_online", x_size = 800, y_size = 600)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -30)

# Create an iso surface.
iso = IsoSurface(scene = s, data_collector = dc)
# Sets the contour number and value.
iso.setValue(contour_number = 1, value = 0.1)
s.render()

