from esys.pyvisi import DataCollector, Map, Scene 
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1000, 
        y_size = 800)

# One DataCollector can be used to specify only one scalar, vector and tensor 
# attribute from a source. If a second scalar, vector or tensor attribute 
# needs to be specified from the same source, a second DataCollector must be 
# created as shown below. In the case below, two scalar attributes 'scalar1' 
# and 'scalar2' were specified.
dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/results.xml")
dc1.setActiveScalar(scalar = "scalar1")

dc2 = DataCollector(source = Source.XML)
dc2.setFileName(file_name = 
        "/home/jongui/results.xml")
dc2.setActiveScalar(scalar = "scalar2")

# Create a map instance for the first viewport.
m1 = Map(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST, 
        lut = Lut.COLOR, outline = True)
m1.setRepresentationToWireframe()

# Create a map instance for the second viewport.
m2 = Map(scene = s, data_collector = dc2, viewport = Viewport.NORTH_WEST, 
        lut = Lut.COLOR, outline = True)
m2.setColor(color = Color.BLUE)

# Create a map instance for the third viewport.
m3 = Map(scene = s, data_collector = dc1, viewport = Viewport.NORTH_EAST, 
        lut = Lut.GREY_SCALE, outline = True)

# Create a map instance the fourth viewport.
m4 = Map(scene = s, data_collector = dc2, viewport = Viewport.SOUTH_EAST, 
        lut = Lut.COLOR, outline = True)
m4.setOpacity(opacity = 0.5)

s.render()
