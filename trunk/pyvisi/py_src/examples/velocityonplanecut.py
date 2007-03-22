from esys.pyvisi import Scene, DataCollector, VelocityOnPlaneCut
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
#dc1.setFileName(file_name = 
dc1.setFileName(file_name = "/home/jongui/data/laurent/slab/source/slab.xml")
#        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")
#dc1.setActiveVector("velocity_cell")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create a velocity instance.
vopc1 = VelocityOnPlaneCut(scene = s, data_collector = dc1, 
        arrow = Arrow.THREE_D, color_mode = ColorMode.VECTOR)
vopc1.setScaleFactor(scale_factor = 100000)
#vopc1.setPlaneToXY(50000)
vopc1.setPlaneToXY()
vopc1.setRatio(3)
vopc1.randomOn()
vopc1.setRatio(25)
#vopc1.setDimension(x = 4, y = 4, z = 4)


s.render()

