from scene import Scene
from datacollector import DataCollector
from camera import Camera
from contour import ContourOnPlane, ContourOnClip
from geo import Transform

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -70)
tf = Transform()
tf.rotateY(angle = 20)

# Creates contours on a clip.
coc = ContourOnClip(scene = s, data_collector = dc, transform = tf)
coc.generateValues(number_contours = 5, min_range = 0.0, max_range = 1.2)
coc.setInsideOutOff()
s.render()

