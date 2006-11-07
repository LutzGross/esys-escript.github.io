from scene import Scene
from datacollector import DataCollector
from camera import Camera
from carpet import Carpet
from geo import Transform

s = Scene(renderer = "vtk_online", x_size = 700, y_size = 600)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../../../data/lutz/seismic/disp.699.vtu")

cam = Camera(scene = s, data_collector = dc)
cam.isometricView()
cam.elevation(angle = -20)
tf = Transform()
tf.xyPlane(offset = -25000)

# Create a carpet for scalar data.
cpt = Carpet(scene = s, data_collector = dc, transform = tf, deform = "Scalar")
# Set the carpet's scale factor.
cpt.setScaleFactor(100000)
s.render()

