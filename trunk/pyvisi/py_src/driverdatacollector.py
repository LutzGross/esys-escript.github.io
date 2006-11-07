from scene import Scene
from datacollector import DataCollector

s = Scene(renderer = "vtk_online", x_size = 500, y_size = 500)
# Create an outline and cube axes for the rendered object. In this case
# the rendered object will be empty as nothing has been specified.
dc = DataCollector(scene = s, outline = True, cube_axes = True)
# The file from which data is to be read.
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")
s.render()
