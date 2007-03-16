from esys.pyvisi import Scene, ImageReader, Image, GlobalPosition,\
DataCollector, Map
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# Create a map instance for the first viewport.
m1 = Map(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST, 
        lut = Lut.COLOR, outline = True)
m1.setOpacity(0.1)

# Create one image reader instance (used in place of data collector).
ir = ImageReader(ImageFormat.JPG)
ir.setImageName(image_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/Flinders_eval.jpg")
	

# Create one image instance.
i = Image(scene = s, image_reader = ir)
#i.setOpacity(opacity = 0.9)
i.translate(0,0,1.)
#i.rotateX(20)
i.setPoint1(GlobalPosition(3,0,0))
i.setPoint2(GlobalPosition(0,3,0))

s.render()

