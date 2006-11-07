from scene import Scene
from datacollector import DataCollector
from camera import Camera
from arrows import ArrowsOnPlane
from geo import Transform

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -50)
tf = Transform()
tf.rotateY(angle = 30)

# Create arrows on a plane.
aop = ArrowsOnPlane(scene = s, data_collector = dc, transform = tf)
aop.setScaleFactor(scale_factor = 0.5)
s.render()

