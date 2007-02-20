from esys.pyvisi import Scene, DataCollector, Ellipsoid, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# Create one ellipsoid instance.
e1 = Ellipsoid(scene = s, data_collector = dc1, tensor = "stress_cell", 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
e1.setScaleFactor(scale_factor = 0.1)
e1.setPhiResolution(10)
e1.setThetaResolution(10)
e1.setDimension(2,2,2)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()
c1.elevation(-20)

s.render()

