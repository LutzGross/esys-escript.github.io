from esys.pyvisi import Scene, DataCollector, Carpet, Camera
from esys.pyvisi.constant import *
import time

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0001.vtk")

cpt = Carpet(scene = s, data_collector = dc1, lut = Lut.COLOR, 
        warp_mode = WarpMode.SCALAR, viewport = Viewport.SOUTH_WEST)
cpt.setPlaneToXY(500000)
cpt.setScaleFactor(0.2)

c = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c.isometricView()

for i in range(1, 100):
    dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.%04d.vtk" % i)

    s.animate()
