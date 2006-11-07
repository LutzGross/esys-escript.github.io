from scene import Scene
from datacollector import DataCollector
from camera import Camera
from map import Map
from colormap import BlueToRed

s = Scene(renderer = "vtk_jpeg", x_size = 1100, y_size = 800)
dc = DataCollector(scene = s, outline = True, cube_axes = False)
# NOTE: Must be executed before the map.
dc.setFileName(file_name = "../../../data/lutz/movie/temp.0.xml")

cam = Camera(scene = s, data_collector = dc)
cam.azimuth(180)
cam.elevation(30)
cam.roll(-20)

br = BlueToRed()
# NOTE: Must be executed before the for loop. Otherwise, results can be
# incorrect.
map = Map(scene = s, data_collector = dc, lut = br)

# Read multiple files. 
for i in range(0 ,80):
	dc.setFileName(file_name = "../../../data/lutz/movie/temp.%d.xml" %i)
	# Save the rendered object from each file into an image.
	#s.saveImage("%02d.jpg" % i)
	# Display the rendered object from each file one after another. 
	s.animate()
