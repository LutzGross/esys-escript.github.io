from scene import Scene
from datacollector import DataCollector
from camera import Camera
from arrows import ArrowsOnPlane, ArrowsOnClip
from geo import Transform

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -60)
tf = Transform()
tf.rotateY(angle = 15)

# Creates arrows on a clip.
aoc = ArrowsOnClip(scene = s, data_collector = dc, transform = tf)
s.render()



