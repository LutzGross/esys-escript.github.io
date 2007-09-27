from esys.pyvisi import Scene, DataCollector, MapOnScalarClip, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = "/home/jongui/data/laurent/slab/source/slab.xml")

# Create a map on scalar clip instance.
mosc1_1 = MapOnScalarClip(scene = s, data_collector = dc1, lut = Lut.GREY_SCALE)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()

s.render()
