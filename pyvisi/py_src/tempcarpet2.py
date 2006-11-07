from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import MapOnPlane
from geo import Transform

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True, cube_axes = True)
dc.setFileName(file_name = "../../../data/lutz/seismic/disp.0.vtu")
#dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")
#dc.setFileName(file_name = "../../../data/laurent/clipper/phi3D.initial.vtk")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(angle = -65)
# Transform is used to change the plane's settings.
tf = Transform()
# Change the plane's default settings.
#tf.rotateX(angle = 15)
#tf.xyPlane(offset = -10000)

# Create a surface map on a plane.
mop = MapOnPlane(scene = s, data_collector = dc, transform = tf)
s.render()
