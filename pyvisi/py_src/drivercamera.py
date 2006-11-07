from scene import Scene
from datacollector import DataCollector
from camera import Camera
from geo import Position

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

# Create a camera. 
cam = Camera(scene = s, data_collector = dc)
# Change the camera's default settings.
cam.setPosition(Position(1.5, 1.5, 10))
cam.setFocalPoint(Position(1, 1, 0.2))
cam.azimuth(angle = 50)

s.render()
