from esys.pyvisi import Scene, DataCollector
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, x_size = 800, y_size = 600, 
        num_viewport = 1)

# Create a data collector intance.
dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")
dc1.setActiveScalar(scalar = "temperature_cell")
dc1.setActiveVector(vector = "velocity_cell")

s.render()
