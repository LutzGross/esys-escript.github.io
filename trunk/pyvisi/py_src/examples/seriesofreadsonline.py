from esys.pyvisi import Scene, DataCollector, MapOnScalarClip, Camera, Map, Contour
from esys.pyvisi.constant import *
import time

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/talus/source/phi_talus_lava.0200.vtu")
dc1.setActiveScalar("lava")

dc2 = DataCollector(source = Source.XML)
dc2.setFileName(file_name = 
        "/home/jongui/data/laurent/talus/source/phi_talus_lava.0200.vtu")

mosc1_1 = Contour(scene = s, data_collector = dc1, lut = Lut.COLOR, 
        viewport = Viewport.SOUTH_WEST)
mosc1_1.generateContours(0)

mosc1_2 = Contour(scene = s, data_collector = dc2, lut = Lut.COLOR, 
        viewport = Viewport.SOUTH_WEST)
mosc1_2.generateContours(0)


# Generate multiple images from multiple files.
for i in range(90, 600):
    dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/talus/source/phi_talus_lava.%04d.vtu" % i)
    dc1.setActiveScalar("lava")
    dc2.setFileName(file_name = 
        "/home/jongui/data/laurent/talus/source/phi_talus_lava.%04d.vtu" % i)

    s.animate()
