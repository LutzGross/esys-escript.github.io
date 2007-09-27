# $Id: arrowPlot3D.py,v 1.1 2005/11/08 08:21:25 paultcochrane Exp $

"""
Example of plotting a 3D vector field
"""

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

#### original vtk code
import vtk

# loading a vtk file as input
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("../../vel-0500.vtk")
reader.Update()

grid = reader.GetOutput()

# grab the model centre and bounds
centre = grid.GetCenter()
bounds = grid.GetBounds()

# grab the norm of the vectors
norm = vtk.vtkVectorNorm()
norm.SetInput(grid)

maxNorm = grid.GetPointData().GetVectors().GetMaxNorm()

# to make arrow glyphs need an arrow source
arrow = vtk.vtkArrowSource()

# the arrows are 3D glyphs so set that up now
glyph = vtk.vtkGlyph3D()
glyph.ScalingOn()
glyph.SetScaleModeToScaleByScalar()
glyph.SetColorModeToColorByScalar()
glyph.SetVectorModeToUseVector()
glyph.SetScaleFactor(0.1/maxNorm)
glyph.SetInput(norm.GetOutput())
glyph.SetSource(arrow.GetOutput())
glyph.ClampingOff()

# set up a stripper to speed up rendering
stripper = vtk.vtkStripper()
stripper.SetInput(glyph.GetOutput())

# make a lookup table for the colour map and invert it (colours look
# better when it's inverted)
lut = vtk.vtkLookupTable()
refLut = vtk.vtkLookupTable()
lut.Build()
refLut.Build()
for j in range(256):
    lut.SetTableValue(j, refLut.GetTableValue(255-j))

# set up the mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInput(stripper.GetOutput())
mapper.SetScalarRange(0,maxNorm)
mapper.SetLookupTable(lut)

# set up the actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# set up the text properties for nice text
font_size = 20
textProp = vtk.vtkTextProperty()
textProp.SetFontSize(font_size)
textProp.SetFontFamilyToArial()
textProp.BoldOff()
textProp.ItalicOff()
textProp.ShadowOff()
textProp.SetColor(0.0, 0.0, 0.0)

# make a title
title = vtk.vtkTextMapper()
title.SetInput("Example 3D arrow/quiver/vector field plot")

# make the title text use the text properties
titleProp = title.GetTextProperty()
titleProp.ShallowCopy(textProp)
titleProp.SetJustificationToCentered()
titleProp.SetVerticalJustificationToTop()
titleProp.BoldOn()

# make the actor for the title
titleActor = vtk.vtkTextActor()
titleActor.SetMapper(title)
titleActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
titleActor.GetPositionCoordinate().SetValue(0.5, 0.95)

# put an outline around the data
outline = vtk.vtkOutlineSource()
outline.SetBounds(bounds)

# make its mapper
outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetInput(outline.GetOutput())

# make its actor
outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(0,0,0)

# set up the renderer and render window
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()

renWin.SetSize(800,600)
renWin.AddRenderer(ren)
ren.SetBackground(1,1,1)

# add the relevant actors
ren.AddActor(actor)
ren.AddActor(titleActor)
ren.AddActor(outlineActor)

cam = ren.GetActiveCamera()
#cam.Azimuth(0)
#cam.Elevation(-90)
cam.Zoom(1.2)
ren.SetActiveCamera(cam)
ren.ResetCameraClippingRange()

# add some axes
axes = vtk.vtkCubeAxesActor2D()
axes.SetInput(grid)
axes.SetCamera(ren.GetActiveCamera())
axes.SetLabelFormat("%6.4g")
axes.SetFlyModeToOuterEdges()
axes.SetFontFactor(0.8)
axes.SetAxisTitleTextProperty(textProp)
axes.SetAxisLabelTextProperty(textProp)
axes.SetXLabel("x")
axes.SetYLabel("y")
axes.SetZLabel("z")
axes.SetNumberOfLabels(5)
axes.GetProperty().SetColor(0,0,0)
ren.AddProp(axes)

# set up stuff for interactive viewing
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

iren.Initialize()
renWin.Render()
iren.Start()

# the WindowToImageFilter is what one uses to save the window to an 
# image file
win2img = vtk.vtkWindowToImageFilter()
win2img.SetInput(renWin)

# set up the PNMWriter as we're saving to png
writer = vtk.vtkPNGWriter()
writer.SetFileName("arrowPlot3D.png")
writer.SetInput(win2img.GetOutput())
writer.Write()

# vim: expandtab shiftwidth=4:
