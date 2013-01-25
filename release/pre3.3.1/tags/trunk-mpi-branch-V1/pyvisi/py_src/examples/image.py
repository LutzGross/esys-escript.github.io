from esys.pyvisi import Scene, ImageReader, Image
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

# Create one image reader instance (used in place of data collector).
ir = ImageReader(ImageFormat.JPG)
ir.setFileName(
        "/home/jongui/trunk/pyvisi/test/python/data_data/Flinders_eval.jpg")

# Create one image instance.
i = Image(scene = s, image_reader = ir)
i.setOpacity(0.7)

s.render()

