from esys.pyvisi import Scene, DataCollector, Map, Light, GlobalPosition, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0271.vtk")

m1 = Map(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()

# Create one light instance.
lig = Light(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
lig.setColor(color = Color.WHITE)
lig.setIntensity(intensity = 1)
lig.setAngle(elevation = -15, azimuth = 0)

# Alternative to using setAngle.
#lig.setPosition(GlobalPosition(0.1,0.1,0.1))
#lig.setFocalPoint(GlobalPosition(0,0,0))

s.render()


