from esys.pyvisi import Scene, DataCollector, Contour, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0271.vtk")

# Create one contour instance.
ctr1 = Contour(scene = s, data_collector = dc1, scalar = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
ctr1.generateContours(1, 0.5, 0.5)

cam1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
cam1.elevation(-40)

s.render()

