from esys.pyvisi import Scene, DataCollector, Velocity 
from esys.pyvisi.constant import *
import time

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0001.vtk")

v = Velocity(scene = s, data_collector = dc1, lut = Lut.COLOR, 
        viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.SCALAR, 
		arrow = Arrow.THREE_D)
v.setScaleFactor(scale_factor = 800000)

for i in range(1, 200):
    print i
    dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.%04d.vtk" % i)

    s.animate()
