from esys.pyvisi import Scene, DataCollector, VelocityOnPlaneClip, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = "/home/jongui/data/laurent/slab/source/slab.xml")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create a velocity on plane clip instance.
vopc1 = VelocityOnPlaneClip(scene = s, data_collector = dc1, 
        arrow = Arrow.THREE_D, color_mode = ColorMode.SCALAR)
vopc1.setScaleFactor(2000000)
vopc1.setPlaneToYZ()
vopc1.rotateY(-70)
vopc1.translate(0,0,0.3)
vopc1.setDimension(7,7,7)

c = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c.isometricView()

s.render()

