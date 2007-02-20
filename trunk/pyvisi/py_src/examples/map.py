from esys.pyvisi import DataCollector, Map, Scene
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0271.vtk")
        #"/home/jongui/data/matt/heat_velocity/source/vel-000659.vtu")

# Create a map instance for the first viewport.
m1 = Map(scene = s, data_collector = dc1, scalar = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
m1.setRepresentationToWireframe()

# Create a map instance for the second viewport.
m2 = Map(scene = s, data_collector = dc1, scalar = None, 
        viewport = Viewport.NORTH_WEST, lut = Lut.COLOR, outline = True)
m2.setColor(Color.BLUE)

# Create a map instance for the third viewport.
m3 = Map(scene = s, data_collector = dc1, scalar = None, 
        viewport = Viewport.NORTH_EAST, lut = Lut.GREY_SCALE, outline = True)

# Create a map instance the fourth viewport.
m4 = Map(scene = s, data_collector = dc1, scalar = None, 
        viewport = Viewport.SOUTH_EAST, lut = Lut.COLOR, outline = True)
m4.setOpacity(0.5)

s.render()
