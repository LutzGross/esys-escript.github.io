from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import Map
from colormap import BlueToRed, RedToBlue

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(-50)
# Create a surface map.
map = Map(scene = s, data_collector = dc)
# Change the LUT to blue-to-red. Default LUT is red-to-blue.
#br = BlueToRed()
#map = Map(scene = s, data_collector = dc, lut = br)

s.render()
