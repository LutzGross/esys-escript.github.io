# $Id: arrowPlot3D.py,v 1.1 2005/11/08 08:22:33 paultcochrane Exp $

"""
Example of plotting a 3D vector field

There are still problems here with getting the angles exactly right....

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


# set up some data to plot
from Numeric import *

dim = 10

# initialise the positions of the vectors
x = zeros((dim,dim), typecode=Float)
y = zeros((dim,dim), typecode=Float)
z = zeros((dim,dim), typecode=Float)

# initialise the vector displacements
# (I may need to rethink how this works in the interface)
dx = zeros((dim,dim), typecode=Float)
dy = zeros((dim,dim), typecode=Float)
dz = zeros((dim,dim), typecode=Float)

# set the positions randomly, and set the displacements to some smaller
# random number but of mean zero instead of distributed between 0 and 1
import random
random.seed()
for i in range(dim):
    for j in range(dim):
        x[i,j] = random.random()
        y[i,j] = random.random()
        z[i,j] = random.random()
        dx[i,j] = (random.random()-0.5)/5.0
        dy[i,j] = (random.random()-0.5)/5.0
        dz[i,j] = (random.random()-0.5)/5.0

#### original povray code

import vtk
import os, sys, re, math
from Numeric import *

# read in the file
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("vel-0500.vtk")
reader.Update()

# get the grid
grid = reader.GetOutput()

# grab the model centre and bounds
centre = grid.GetCenter()
bounds = grid.GetBounds()

# try and extract the vector norm
norm = vtk.vtkVectorNorm()
norm.SetInput(grid)

maxNorm = grid.GetPointData().GetVectors().GetMaxNorm()

### extract the relevant grid data

# the points
points = grid.GetPoints()
numPoints = points.GetNumberOfPoints()
x = zeros(numPoints, typecode=Float)
y = zeros(numPoints, typecode=Float)
z = zeros(numPoints, typecode=Float)
for i in range(numPoints):
    x[i], y[i], z[i] = points.GetPoint(i)

# the data at the points
data = grid.GetPointData().GetVectors()
vx = zeros(numPoints, typecode=Float)
vy = zeros(numPoints, typecode=Float)
vz = zeros(numPoints, typecode=Float)
vNorm = zeros(numPoints, typecode=Float)
for i in range(numPoints):
    vx[i], vy[i], vz[i] = data.GetTuple3(i)
    vNorm[i] = math.sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])

# make a lookup table for the colour map and invert it (colours look
# better when it's inverted)
lut = vtk.vtkLookupTable()
refLut = vtk.vtkLookupTable()
lut.Build()
refLut.Build()
for j in range(256):
    lut.SetTableValue(j, refLut.GetTableValue(255-j))

# get the colours
r = zeros(numPoints, typecode=Float)
g = zeros(numPoints, typecode=Float)
b = zeros(numPoints, typecode=Float)
for i in range(numPoints):
    r[i], g[i], b[i] = lut.GetColor(vNorm[i]/maxNorm)

### generate the pov file

pov = open("arrowPlot3D.pov", "w")

pov.write("#include \"colors.inc\"\n")
pov.write("#include \"shapes.inc\"\n")
pov.write("#include \"textures.inc\"\n")

pov.write("camera {\n")
pov.write("  location <%f, %f, -2.5>\n" % (centre[0], centre[1]))
pov.write("  look_at <%f, %f, -%f>\n" % (centre[0], centre[1], centre[2]))
pov.write("}\n")

pov.write("light_source {\n")
pov.write("  <0, 0, -3>\n")
pov.write("  colour White\n")
pov.write("}\n")

pov.write("#declare Arrow = union {\n")
pov.write("  cone {\n")
pov.write("    <0, 0, 0>, 0.3\n")
pov.write("    <1, 0, 0>, 0.0\n")
pov.write("  }\n")
    
pov.write("  cylinder {\n")
pov.write("    <-1, 0, 0>\n")
pov.write("    <0, 0, 0>,\n")
pov.write("    0.15\n")
pov.write("  }\n")
pov.write("}\n")

for i in range(numPoints):
    pov.write("object {\n")
    scale = 0.05*vNorm[i]/maxNorm
    if scale < 1e-8:
	scale = 1e-7
    pov.write("  Arrow scale %g " % scale)
    pov.write("rotate <%f, %f, %f> " % (vx[i], vy[i], vz[i]))
    pov.write("translate <%f, %f, -%f> " % (x[i], y[i], z[i]))
    pov.write("pigment { colour <%f, %f, %f> }\n" % (r[i], g[i], b[i]))
    pov.write("}\n")

pov.close()

### generate the ini file

# open the ini file to write to
ini = open("arrowPlot3D.ini", "w")

# the output resolution
ini.write("Width=640\n")
ini.write("Height=480\n")

# anti-aliasing settings
ini.write("Antialias=on\n")

# generate png files
ini.write("Output_File_Type=N\n")

# the name of the input pov file
ini.write("Input_File_Name=arrowPlot3D.pov\n")

# pause when done
ini.write("Pause_When_Done=on\n")

# close the file
ini.close()

# run povray on the file
result = os.system("povray arrowPlot3D.ini")
if result != 0:
    raise SystemError, "Povray execution failed"
else:
    # clean up a bit
    os.unlink("arrowPlot3D.pov")
    os.unlink("arrowPlot3D.ini")

# vim: expandtab shiftwidth=4:
