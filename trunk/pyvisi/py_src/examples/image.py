from esys.pyvisi import Scene, ImageReader, Image, GlobalPosition
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

# Create one image reader instance (used in place of data collector).
ir = ImageReader(ImageFormat.JPG)
ir.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/Flinders_eval.jpg")

# Create one image instance.
i = Image(scene = s, image_reader = ir)
i.setOpacity(opacity = 0.9)
i.setPosition(GlobalPosition(-600,50.9,0.5))

s.render()

