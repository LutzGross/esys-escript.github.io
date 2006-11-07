from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import MapOnPlane
from geo import Transform
import time

s = Scene(renderer = "vtk_jpeg", x_size = 1000, y_size = 800)

dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.isometricView()

# Demonstrates a simple animation with two planes on a map. 
tf = Transform()
tf.xzPlane()

tf2 = Transform()
tf2.yzPlane()

mop = MapOnPlane(scene = s, data_collector = dc, transform = tf)
mop2 = MapOnPlane(scene = s, data_collector = dc, transform = tf2)

for i in range(0,60):
	tf.translate(0,0.05,0)
	tf2.translate(0.05,0,0)
	# This creates on the fly animation. No window interaction can occur.
	s.animate()
	#s.saveImage("%02d.jpg" % i)
	time.sleep(0.05)

