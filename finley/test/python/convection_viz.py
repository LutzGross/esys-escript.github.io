
# Create an MPEG movie called movie.mpg from the output of convection.py using pyvisi
# Usage: python convection_viz.py

# View movie with: mplayer movie.mpg
# Create an animated GIF "out.gif" with: mplayer -vo gif89a movie.mpg

# The output of convection.py consists of VTK files state.i1.vtu ... state.i2.vtu
# Specify range of file names:   state.i1.vtu ... state.i2.vtu
i1 = 0
i2 = 271


################################################################################
# You shouldn't have to change anything below here
################################################################################


from esys.pyvisi import Camera, Contour, DataCollector, Legend, LocalPosition, Map, Movie, Scene, Velocity, VelocityOnPlaneCut
from esys.pyvisi.constant import *
import os, sys

# List the input files (which are in .vtu format)
filList = []
for i in range(i1, i2+1):
  filList.append("state.%d.vtu" % (i))

print "Making a movie using files:  ", filList[0], "...", filList[len(filList)-1], "\n"

# Make sure all input files exist before starting
for filName in filList:
  # Check that the name is a valid file
  if not os.path.isfile(filName):
    print "%s: Invalid input file '%s'\n" % (sys.argv[0], filName)
    sys.exit(1)

# Now create the .jpg images for each input VTK file
imgList = []
for filName in filList:
  imgName = filName + ".jpg"

  # Check that the name is a valid file
  if not os.path.isfile(filName):
    print "%s: Input file '%s' has disappeared\n" % (sys.argv[0], filName)
    sys.exit(1)

  print "Reading %s and writing %s" % (filName, imgName)
  imgList.append(imgName)

  # Create a pyvisi Scene
  s = Scene(renderer = Renderer.OFFLINE_JPG, num_viewport = 1, x_size = 400, y_size = 400)

  # Create a DataCollector reading a XML file
  dc1 = DataCollector(source = Source.XML)
  dc1.setFileName(file_name = filName)
  dc1.setActiveScalar(scalar = "T")
  dc1.setActiveVector(vector = "v")

  # Create a Contour
  ctr1 = Contour(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, cell_to_point = False, outline = True)
  # ctr1.generateContours(contours = 1, lower_range=0.5, upper_range=0.5) # lower_range controls top isosurface in the case of convection.py
  ctr1.generateContours(contours = 2, lower_range=0.3, upper_range=0.9) # lower_range controls top isosurface in the case of convection.py

  # Create VelocityOnPlaneCut at 0.9 to show the direction of movement where the continents would be
  if False:
    vopc1 = VelocityOnPlaneCut(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST, color_mode = ColorMode.VECTOR, arrow = Arrow.THREE_D, lut = Lut.COLOR, cell_to_point = False, outline = True)
    vopc1.setScaleFactor(scale_factor = 70.0)	# Set the scale factor for the length of the arrows
    vopc1.setPlaneToXY(offset = 0.9)		# Set the height of the plane orthogonal to the z-axis
    vopc1.setRatio(2)				# Mask every Nth point
    vopc1.randomOn()

  # Create a Camera
  cam1 = Camera(scene = s, viewport = Viewport.SOUTH_WEST)
  cam1.elevation(angle = -60)

  # Render the jpg image for this input file
  s.render(image_name = imgName)


# Finally, create the movie from the .jpg files
mov = Movie()
mov.imageList(input_directory = ".", image_list = imgList)
mov.makeMovie("movie.mpg")

print "\nCreated movie.mpg\n"

