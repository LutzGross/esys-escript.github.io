from esys.pyvisi import DataCollector, Velocity, Scene
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# Create a velocity instance in the first viewport.
v1 = Velocity(scene = s, data_collector = dc1, vector = None, 
        viewport = Viewport.SOUTH_WEST, arrow = Arrow.THREE_D, 
        color_mode = ColorMode.VECTOR, lut = Lut.COLOR, outline = True)
v1.setRepresentationToWireframe()
v1.setScaleFactor(0.2)
v1.setScaleModeByScalar()
v1.setDimension(2,2,2)

# Create a velocity instance in the second viewport.
v2 = Velocity(scene = s, data_collector = dc1, vector = None, 
        viewport = Viewport.NORTH_WEST, arrow = Arrow.THREE_D, 
        color_mode = ColorMode.SCALAR, lut = Lut.COLOR, outline = True)
v2.setColor(Color.BLUE)
v2.setScaleFactor(0.2)

# Create a velocity instance in the third viewport.
v3 = Velocity(scene = s, data_collector = dc1, vector = None, 
        viewport = Viewport.NORTH_EAST, arrow = Arrow.TWO_D, 
        color_mode = ColorMode.VECTOR, lut = Lut.GREY_SCALE, outline = True)
v3.setScaleFactor(0.2)
v3.setDimension(4,4,4)

# Create a velocity instance in the fourth viewport.
v4 = Velocity(scene = s, data_collector = dc1, vector = None, 
        viewport = Viewport.SOUTH_EAST,  arrow = Arrow.TWO_D, 
        color_mode = ColorMode.SCALAR, lut = Lut.COLOR, outline = True)
v4.setOpacity(0.5)
v4.setScaleFactor(0.2)

s.render()
