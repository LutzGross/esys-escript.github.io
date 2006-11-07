from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import MapOnPlane, MapOnClip
from geo import Transform

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.roll(angle = -90)
cam.elevation(angle = -65)
tf = Transform()
tf.rotateX(15)

# Create a map on a clip.
moc = MapOnClip(scene = s, data_collector = dc, transform = tf)
# Clip the other side of the rendered object.
#moc.setInsideOutOff()

s.render()


