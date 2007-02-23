from esys.pyvisi import Scene, DataCollector, MapOnPlaneClip, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# NOTE: There is a difference between performing rotation then followed by 
# translation and performing translation then followed by rotation.

# Create three map on clip instances.
mopc1_1 = MapOnPlaneClip(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST)
mopc1_1.setPlaneToXY()
mopc1_1.rotateX(angle = 5)

mopc1_2 = MapOnPlaneClip(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST)
mopc1_2.setPlaneToYZ(offset = 2.5)
mopc1_2.setOpacity(opacity = 0.5)

mopc1_3 = MapOnPlaneClip(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST)
mopc1_3.setPlaneToXZ()
mopc1_3.rotateX(angle = -40)
mopc1_3.translate(x_offset = 0, y_offset = 2.2, z_offset = 0)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()

s.render()
