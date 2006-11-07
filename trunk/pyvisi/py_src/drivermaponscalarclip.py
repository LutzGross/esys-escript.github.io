from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import MapOnClip, MapOnScalarClip

s = Scene(renderer = "vtk_online", x_size = 500, y_size = 600)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../../../data/laurent/clipper/phi3D.initial.vtk")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -60)
cam.azimuth(angle = -20)

# Create a map on a scalar clip.
moc = MapOnScalarClip(scene = s, data_collector = dc)
s.render()


