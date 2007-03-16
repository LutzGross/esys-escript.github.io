from esys.pyvisi import Scene, DataCollector, Map 
from esys.pyvisi.constant import *
import time

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/doc/examples/T.1.xml")

m = Map(scene = s, data_collector = dc1, lut = Lut.COLOR, 
        viewport = Viewport.SOUTH_WEST)

# Generate multiple images from multiple files.
for i in range(1, 50):
    dc1.setFileName(file_name = 
        "/home/jongui/trunk/doc/examples/T.%d.xml" % i)

    s.animate()
