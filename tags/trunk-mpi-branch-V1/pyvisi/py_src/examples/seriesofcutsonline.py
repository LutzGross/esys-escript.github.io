from esys.pyvisi import Scene, DataCollector, MapOnPlaneCut, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

# Create two data collector instances.
dc1 = DataCollector(source = Source.XML)

dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# NOTE: There is a difference between performing rotation then followed by 
# translation and performing translation then followed by rotation.

# Create a map on plane cut instance for the first viewport.
mopc1 = MapOnPlaneCut(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST)
mopc1.setPlaneToYZ(0.1)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()

# Animate multiple cuts from the translation.
for i in range(0, 60):
    s.animate()
    mopc1.translate(0.05,0,0)
