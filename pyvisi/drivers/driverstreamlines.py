from scene import Scene
from datacollector import DataCollector
from camera import Camera
from streamlines import StreamLines 

s = Scene(renderer = "vtk_online", x_size = 1000, y_size = 800)
dc = DataCollector(scene = s, outline = True)
dc.setFileName(file_name = "../test/python/data_data/interior_3D.xml")

cam = Camera(scene = s, data_collector = dc)
cam.elevation(-50)

# Create streamlines
sl = StreamLines(scene = s, data_collector = dc)
# Change the streamlines default accuracy.
sl.setAccuracy(0.2)
s.render()

