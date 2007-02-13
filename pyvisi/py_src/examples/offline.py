from esys.pyvisi import Scene, DataCollector, Map
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.OFFLINE_JPG, x_size = 800, y_size = 600)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

m1 = Map(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)

# Save the rendered object into an image.
s.saveImage("offline_image.jpg")

