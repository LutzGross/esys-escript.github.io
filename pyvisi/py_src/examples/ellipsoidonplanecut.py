from esys.pyvisi import Scene, DataCollector, EllipsoidOnPlaneCut, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create one ellipsoid on plance cut instance.
eopc1 = EllipsoidOnPlaneCut(scene = s, data_collector = dc1, tensor = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
eopc1.setScaleFactor(scale_factor = 0.1)
eopc1.setPlaneToXY()
eopc1.rotateX(angle = 20)
eopc1.setDimension(x = 2, y = 2, z = 2)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()
c1.elevation(angle = -20)

s.render()

