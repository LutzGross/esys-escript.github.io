from esys.pyvisi import Scene, DataCollector, StreamLine, Camera, GlobalPosition
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# Create one streamline instance for the first viewport.
sl1 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True, 
        color_mode = ColorMode.SCALAR)
sl1.setTubeRadius(0.01)

cam1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
cam1.elevation(-40)

# Create one streamline instance for the second viewport.
sl2 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.NORTH_WEST, lut = Lut.COLOR, outline = True)
sl2.setPointSourceRadius(0.4)
sl2.setPointSourceNumberOfPoints(5)
sl2.setTubeRadiusToVaryByScalar()

cam2 = Camera(scene = s, data_collector = dc1, viewport = Viewport.NORTH_WEST)
cam2.elevation(-40)

# Create one streamline instance for the third viewport.
sl3 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.NORTH_EAST, lut = Lut.COLOR, outline = True, 
        color_mode = ColorMode.SCALAR)
sl3.setPointSourceCenter(GlobalPosition(0.7,0.7,0.2))

cam3 = Camera(scene = s, data_collector = dc1, viewport = Viewport.NORTH_EAST)
cam3.elevation(-60)

# Create one streamline instance for the fourth viewport. 
sl4 = StreamLine(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_EAST, lut = Lut.COLOR, outline = True)
sl4.setPointSourceCenter(GlobalPosition(0.7,0.7,0.2))
sl4.setMaximumPropagationTime(100)
sl4.setStepLength(0.5)
sl4.setIntegrationStepLength(2.0)

cam4 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_EAST)
cam4.elevation(-60)

s.render()

