from esys.pyvisi import Scene, DataCollector, ContourOnPlaneCut, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0271.vtk")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create one contour on plane cut instance.
ctropc1 = ContourOnPlaneCut(scene = s, data_collector = dc1, scalar = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
ctropc1.generateContours(contours = 8)
ctropc1.setPlaneToXY(offset = 200000)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.elevation(angle = -45)

s.render()

