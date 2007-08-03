from esys.pyvisi import Scene, DataCollector, EllipsoidOnPlaneClip, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create on ellipsoid on plane clip instance.
eopc = EllipsoidOnPlaneClip(scene = s, data_collector = dc1, tensor = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
eopc.setScaleFactor(0.1)
eopc.setPlaneToXY()
eopc.rotateX(20)
eopc.setDimension(2,2,2)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.bottomView()
c1.azimuth(-90)
c1.elevation(20)

s.render()

