from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import Map
from light import Light
from constants import *

s = Scene(renderer = "vtk_online", x_size = 700, y_size = 600)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.isometricView()
# Creates a light source.
lig = Light(scene = s, data_collector = dc)
lig.setColor(PURPLE)

map = Map(scene = s, data_collector = dc)
s.render()
