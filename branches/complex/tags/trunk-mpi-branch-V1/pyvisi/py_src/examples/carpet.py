from esys.pyvisi import Scene, DataCollector, Carpet, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# Create one carpet instance.
cpt1 = Carpet(scene = s, data_collector = dc1, warp_mode = WarpMode.SCALAR, 
        lut = Lut.COLOR)
cpt1.setPlaneToXY(0.5)
cpt1.setScaleFactor(0.8)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()

s.render()

