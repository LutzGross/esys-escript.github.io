from esys.pyvisi import Scene, DataCollector, StreamLine, Camera, GlobalPosition
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/results.xml")

# Create one streamline instance for the first viewport.
sl1 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True, 
        color_mode = ColorMode.SCALAR, scalar = "scalar2")
sl1.setTubeRadius(radius = 0.01)

cam1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
cam1.elevation(angle = -40)

# Create one streamline instance for the second viewport.
sl2 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.NORTH_WEST, lut = Lut.COLOR, outline = True)
sl2.setPointSourceRadius(radius = 0.4)
sl2.setPointSourceNumberOfPoints(points = 5)
sl2.setTubeRadiusToVaryByScalar()

cam2 = Camera(scene = s, data_collector = dc1, viewport = Viewport.NORTH_WEST)
cam2.elevation(angle = -40)

# Create one streamline instance for the third viewport.
sl3 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.NORTH_EAST, lut = Lut.COLOR, outline = True, 
        color_mode = ColorMode.SCALAR, scalar = "scalar1")
sl3.setPointSourceCenter(GlobalPosition(x_coor = 0.2, y_coor = 0.2, 
        z_coor = 0.2))

cam3 = Camera(scene = s, data_collector = dc1, viewport = Viewport.NORTH_EAST)
cam3.elevation(angle = -60)

# Create one streamline instance for the fourth viewport. 
sl4 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_EAST, lut = Lut.COLOR, outline = True)
sl4.setPointSourceCenter(GlobalPosition(x_coor = 0.7, y_coor = 0.7, 
        z_coor = 0.2))
sl4.setMaximumPropagationTime(100)

cam4 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_EAST)
cam4.elevation(angle = -60)

s.render()

