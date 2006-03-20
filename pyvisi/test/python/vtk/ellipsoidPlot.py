# $Id: ellipsoidPlot.py,v 1.1 2005/09/21 00:37:33 paultcochrane Exp $

"""
Example of plotting ellipsoids (useful for visualising tensors) with pyvisi 
"""

import vtk

# reverse the order of the colourmap (looks better)
lut = vtk.vtkLookupTable()
refLut = vtk.vtkLookupTable()
lut.Build()
refLut.Build()
for i in range(256):
    lut.SetTableValue(i, refLut.GetTableValue(255-i))

# set up the renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(640,480)
ren.SetBackground(1,1,1)  # white

# load a vtk file as input
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("stress22.vtk")
reader.Update()

# grab the grid of the data
grid = reader.GetOutput()

# convert the cell data to point data
c2p = vtk.vtkCellDataToPointData()
c2p.SetInput(grid)

# now extract the tensor components
extract = vtk.vtkExtractTensorComponents()
extract.SetInput(c2p.GetOutput())
extract.SetScalarModeToEffectiveStress()
extract.ExtractScalarsOn()
extract.PassTensorsToOutputOn()
extract.ScalarIsEffectiveStress()

extractGrid = extract.GetOutput()
extractGrid.Update()
extractScalarRange = extractGrid.GetPointData().GetScalars().GetRange()

# make a sphere source for the glyphs
sphere = vtk.vtkSphereSource()
sphere.SetThetaResolution(6)
sphere.SetPhiResolution(6)
sphere.SetRadius(0.5)

# make tensor glyphs
glyph = vtk.vtkTensorGlyph()
glyph.SetSource(sphere.GetOutput())
glyph.SetInput(extractGrid)
glyph.SetColorModeToScalars()
glyph.ScalingOn()
glyph.SetMaxScaleFactor(5.0)
glyph.SetScaleFactor(1.0)
glyph.ClampScalingOn()

# make a stripper for faster rendering
stripper = vtk.vtkStripper()
stripper.SetInput(glyph.GetOutput())

# make the normals of the data
normals = vtk.vtkPolyDataNormals()
normals.SetInput(stripper.GetOutput())

# make the mapper for the data
mapper = vtk.vtkPolyDataMapper()
mapper.SetInput(normals.GetOutput())
mapper.SetLookupTable(lut)
mapper.SetScalarRange(extractScalarRange)

# make the actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# add the actor to be rendered
ren.AddActor(actor)

# set up text properties
textProp = vtk.vtkTextProperty()
textProp.SetFontFamilyToArial()
textProp.BoldOff()
textProp.ItalicOff()
textProp.ShadowOff()
textProp.SetColor(0.0, 0.0, 0.0)

# make a title
title = vtk.vtkTextMapper()
title.SetInput("Example ellipsoid plot")

# make the title text use the text properties
titleProp = title.GetTextProperty()
titleProp.ShallowCopy(textProp)
titleProp.SetJustificationToCentered()
titleProp.SetVerticalJustificationToTop()
titleProp.SetFontSize(20)
titleProp.BoldOn()

# make the actor for the title
titleActor = vtk.vtkTextActor()
titleActor.SetMapper(title)
titleActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
titleActor.GetPositionCoordinate().SetValue(0.5, 0.95)

ren.AddActor(titleActor)

# add a scalar bar
scalarBar = vtk.vtkScalarBarActor()
scalarBar.SetLookupTable(lut)
scalarBar.SetWidth(0.1)
scalarBar.SetHeight(0.8)
scalarBar.SetPosition(0.9, 0.15)

# set up the label text properties 
scalarBarTextProp = scalarBar.GetLabelTextProperty()
scalarBarTextProp.ShallowCopy(textProp)
scalarBarTextProp.SetFontSize(10)

ren.AddActor(scalarBar)

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
writer.SetFileName("ellipsoidPlot.png")
writer.SetInput(win2img.GetOutput())
writer.Write()

# vim: expandtab shiftwidth=4:

