from esys.pyvisi import DataCollector, Map, Scene, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1000, 
        y_size = 800)

dc = DataCollector(source = Source.XML)
dc.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0271.vtk")

m = Map(scene = s, data_collector = dc, scalar = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
m.setRepresentationToWireframe()

# Create a camera instance for the first viewport.
c = Camera(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST)
c.azimuth(angle = 50)

m2 = Map(scene = s, data_collector = dc, scalar = None, 
        viewport = Viewport.NORTH_WEST, lut = Lut.COLOR, outline = True)
m2.setColor(color = Color.BLUE)

# Create a camera instance for the second viewport.
c2 = Camera(scene = s, data_collector = dc, viewport = Viewport.NORTH_WEST)
c2.elevation(angle = 50)

m3 = Map(scene = s, data_collector = dc, scalar = None, 
        viewport = Viewport.NORTH_EAST, lut = Lut.GREY_SCALE, outline = True)

# Create a camera instance for the third viewport.
c3 = Camera(scene = s, data_collector = dc, viewport = Viewport.NORTH_EAST)
c3.isometricView()

m4 = Map(scene = s, data_collector = dc, scalar = None, 
        viewport = Viewport.SOUTH_EAST, lut = Lut.COLOR, outline = True)
m4.setOpacity(opacity = 0.7)

# Create a camera instance for the fourth viewport.
c4 = Camera(scene = s, data_collector = dc, viewport = Viewport.SOUTH_EAST)
c4.leftView()

s.render()
