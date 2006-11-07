from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import MapOnPlane
from geo import Transform

s = Scene(renderer = "vtk_jpeg", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -55)

tf = []
for i in range(0,3):
	tf.append(Transform())
	eval("tf[%d].xzPlane(offset = %d)" % (i,i))
	eval("MapOnPlane(scene = s, data_collector = dc, transform = tf[%d])" % i)

s.render()
# Perform offline rendering  and saving the rendered object as an image.
#s.saveImage("multipleplanes.jpg")


