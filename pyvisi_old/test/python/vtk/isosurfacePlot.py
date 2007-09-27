# $Id: isosurfacePlot.py,v 1.1 2005/11/08 06:08:17 paultcochrane Exp $

"""
Example of plotting a set of isosurfaces
"""

# set up some data to plot
from Numeric import *

import vtk

# load a vtk file as input
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("../../temp-0500.vtk")
reader.Update()

# need to do a delaunay 3D here to get decent looking isosurfaces
del3D = vtk.vtkDelaunay3D()
del3D.SetInput(reader.GetOutput())
del3D.SetOffset(2.5)
del3D.SetTolerance(0.001)
del3D.SetAlpha(0.0)

# set up a contour filter
cont = vtk.vtkContourGrid()
cont.SetInput(del3D.GetOutput())
cont.GenerateValues(5, 0.25, 0.75)  # need to automate this!!
cont.ComputeScalarsOn()

# get the model centre and bounds
centre = reader.GetOutput().GetCenter()
bounds = reader.GetOutput().GetBounds()

# set up the mapper
mapper = vtk.vtkDataSetMapper()
mapper.SetInput(cont.GetOutput())
mapper.ScalarVisibilityOn()

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
title.SetInput("Example isosurface plot")

# make the title text use the text properties
titleProp = title.GetTextProperty()
titleProp.ShallowCopy(textProp)
titleProp.SetJustificationToCentered()
titleProp.SetVerticalJustificationToTop()

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

# make a lookup table for the colour map and invert it (colours look
# better when it's inverted)
lut = vtk.vtkLookupTable()
refLut = vtk.vtkLookupTable()
lut.Build()
refLut.Build()
for j in range(256):
    lut.SetTableValue(j, refLut.GetTableValue(255-j))

# set up the renderer and render window
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()

renWin.SetSize(800,600)
renWin.AddRenderer(ren)
ren.SetBackground(1,1,1)

# add some axes
axes = vtk.vtkCubeAxesActor2D()
axes.SetInput(reader.GetOutput())
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

# add the relevant actors
ren.AddActor(actor)
ren.AddActor(titleActor)
ren.AddActor(outlineActor)
ren.AddActor(axes)

### need to make a camera object in pyvisi!!!
ren.ResetCamera()
cam = ren.GetActiveCamera()
cam.Azimuth(0)
cam.Elevation(-90)
cam.SetFocalPoint(centre)
ren.SetActiveCamera(cam)
ren.ResetCameraClippingRange()

# play around with lighting
light = vtk.vtkLight()
light.SetFocalPoint(centre)
light.SetPosition(centre[0]-bounds[1],
	centre[1]-bounds[3],
	centre[2]+bounds[5])
ren.AddLight(light)

light2 = vtk.vtkLight()
light2.SetFocalPoint(centre)
light2.SetPosition(centre[0]+bounds[1],
	 centre[1]+bounds[3],
	 centre[2]-bounds[5])
ren.AddLight(light2)

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

# set up the PNGWriter as we're saving to png
writer = vtk.vtkPNGWriter()
writer.SetFileName("isosurfacePlot.png")
writer.SetInput(win2img.GetOutput())
writer.Write()

# vim: expandtab shiftwidth=4:
