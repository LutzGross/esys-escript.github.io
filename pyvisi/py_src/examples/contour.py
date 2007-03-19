from esys.pyvisi import Scene, DataCollector, Contour, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")
dc1.setActiveScalar(scalar = "temperature_cell")

# Create one contour instance.
ctr1 = Contour(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST,
        lut = Lut.COLOR, outline = True)
#ctr1.generateContours(contours = 1, lower_range = 0.5, upper_range = 0.5)

cam1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
cam1.elevation(angle = -40)

s.render()

