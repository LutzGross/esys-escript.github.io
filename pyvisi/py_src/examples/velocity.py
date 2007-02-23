from esys.pyvisi import DataCollector, Velocity, Scene
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/results.xml")

dc2 = DataCollector(source = Source.XML)
dc2.setFileName(file_name = 
        "/home/jongui/results.xml")

dc3 = DataCollector(source = Source.XML)
dc3.setFileName(file_name = 
        "/home/jongui/results.xml")

# Create a velocity instance in the first viewport.
v1 = Velocity(scene = s, data_collector = dc1, vector = "vector", 
        scalar = "scalar2", viewport = Viewport.SOUTH_WEST, 
        arrow = Arrow.THREE_D, color_mode = ColorMode.VECTOR, 
        lut = Lut.COLOR, outline = True)
v1.setRepresentationToWireframe()
v1.setScaleFactor(scale_factor = 0.3)
v1.setScaleModeByScalar()
v1.setDimension(x = 2, y = 2, z = 2)

# Create a velocity instance in the second viewport.
v2 = Velocity(scene = s, data_collector = dc2, vector = "vector", 
        scalar = "scalar1", viewport = Viewport.NORTH_WEST, 
        arrow = Arrow.THREE_D, color_mode = ColorMode.VECTOR, 
        lut = Lut.COLOR, outline = True)
v2.setColor(color = Color.BLUE)
v2.setScaleModeByScalar()
v2.setScaleFactor(scale_factor = 0.3)

# Create a velocity instance in the third viewport.
v3 = Velocity(scene = s, data_collector = dc2, vector = "vector", 
        scalar = "scalar1", viewport = Viewport.NORTH_EAST, 
        arrow = Arrow.TWO_D, color_mode = ColorMode.SCALAR, 
        lut = Lut.COLOR, outline = True)
v3.setScaleFactor(scale_factor = 0.2)
v3.setScaleModeByVector()
v3.setDimension(x = 1, y = 1, z = 1)

# Create a velocity instance in the fourth viewport.
v4 = Velocity(scene = s, data_collector = dc3, vector = "vector2", 
        scalar = "scalar1", viewport = Viewport.SOUTH_EAST,  
        arrow = Arrow.TWO_D, color_mode = ColorMode.SCALAR, 
        lut = Lut.COLOR, outline = True)
v4.setOpacity(opacity = 0.5)
v4.setScaleModeByVector()
v4.setScaleFactor(scale_factor = 0.2)

s.render()
