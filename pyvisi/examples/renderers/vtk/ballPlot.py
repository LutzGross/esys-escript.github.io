# $Id: ballPlot.py,v 1.2 2005/05/25 05:39:19 paultcochrane Exp $

"""
Example of plotting spheres with pyvisi 
"""

import os, sys
import random

# set up some data to plot
from Numeric import *

# the three axes in space
# this will give us 10 particles (_not_ 1000)
x = arange(10, typecode=Float)
y = arange(10, typecode=Float)
z = arange(10, typecode=Float)

# 3D position information
posArray = []
for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            posArray.append( (x[i], y[j], z[k]) )

# radius information
random.seed()
radiiArray = zeros(len(x)*len(y)*len(z), typecode=Float)
for i in range(len(x)*len(y)*len(z)):
    radiiArray[i] = random.random()*0.8

# tag information
random.seed()
tagsArray = zeros(len(x)*len(y)*len(z), typecode=Int)
for i in range(len(x)*len(y)*len(z)):
    tagsArray[i] = int(random.random()*10)

#### original vtk code

import vtk
import os

    # we're loading a file; make sure it exists first
if os.path.exists("../../cp_test_0.xml"):
    pass
else:
    raise IOError, "File not found"

# create the reader of the file
_reader = vtk.vtkXMLUnstructuredGridReader()
_reader.SetFileName("../../cp_test_0.xml")
_reader.Update()

# read the output into an unstructured grid
_grid = _reader.GetOutput()

# grab the radius data for the radii of the balls
_radii = _grid.GetPointData().GetScalars("radius")

# grab the tag data and use for colouring the balls
_tags = _grid.GetPointData().GetScalars("particleTag")

# work out dynamically the number of different tags so that can use this
# information to automatically set the scalar range for colouring
_numPoints = _tags.GetNumberOfTuples()
_valueDict = {}
for i in range(_numPoints):
    _tagValue = _tags.GetValue(i)
    _valueDict[_tagValue] = 1

_numTags = len(_valueDict.keys())

_tagValues = _valueDict.keys()
_tagValues.sort()

# count the number of tags, and make an evenly spaced array of points
# between zero and one, then use these as the scalars to colour by
_scaledTags = vtk.vtkFloatArray()
_scaledTags.SetNumberOfTuples(_numPoints)
_scaledTags.SetNumberOfComponents(1)
_scaledTags.SetName("tags")
for i in range(_numPoints):
    _tagValue = _tags.GetValue(i)
    for j in range(_numTags):
        if _tagValues[j] == _tagValue:
            _scaledTags.InsertTuple1(i, float(j)/float(_numTags-1))

# in vtk 4.2 have to set up an array of two components to get
# the data through the glyph object to the mapper so create this now
_data = vtk.vtkFloatArray()
_data.SetNumberOfComponents(3)
_data.SetNumberOfTuples(_radii.GetNumberOfTuples())
_data.CopyComponent(0, _radii, 0)
_data.CopyComponent(1, _tags, 0)
_data.CopyComponent(2, _scaledTags, 0)
_data.SetName("data")

# add the data array to the grid (again, only specific to vtk 4.2 and 4.4)
_grid.GetPointData().AddArray(_data)

# make the data the active scalars
_grid.GetPointData().SetActiveScalars("data")

# to make sphere glyphs need a sphere source
_sphere = vtk.vtkSphereSource()
_sphere.SetRadius(1.0)  # set to 0.5 by default in vtk, we need it to be 1.0
_sphere.SetThetaResolution(5)  # chunky, but will be good in large models
_sphere.SetPhiResolution(5)

# the spheres are 3D glyphs so set that up now
_glyph = vtk.vtkGlyph3D()
_glyph.ScalingOn()
_glyph.SetScaleModeToScaleByScalar() # scale by the radii scalars
_glyph.SetColorModeToColorByScalar() # colour by the tag scalars
_glyph.SetScaleFactor(1.0)  # just in case
_glyph.SetInput(_grid)
_glyph.SetSource(_sphere.GetOutput())  # grab the 3D glyph to use
_glyph.ClampingOff()  # so we even see the tiny spheres

# set up a stripper (this will speed up rendering a lot)
_stripper = vtk.vtkStripper()
_stripper.SetInput(_glyph.GetOutput())

# set up the mapper
_mapper = vtk.vtkPolyDataMapper()
_mapper.SetInput(_stripper.GetOutput())
_mapper.ScalarVisibilityOn()
_mapper.ColorByArrayComponent("data", 2)
_mapper.SetScalarRange(0, 1)  

# set up the actor
_actor = vtk.vtkActor()
_actor.SetMapper(_mapper)

# text properties
_font_size = 14
_textProp = vtk.vtkTextProperty()
_textProp.SetFontSize(_font_size)
_textProp.SetFontFamilyToArial()
_textProp.BoldOff()
_textProp.ItalicOff()
_textProp.ShadowOff()

# add a title
_titleMapper = vtk.vtkTextMapper()
_title = "Example ball plot"
_titleMapper.SetInput(_title)

_titleProp = _titleMapper.GetTextProperty()
_titleProp.ShallowCopy(_textProp)
_titleProp.SetJustificationToCentered()
_titleProp.SetVerticalJustificationToTop()

# set up the text actor
_titleActor = vtk.vtkTextActor()
_titleActor.SetMapper(_titleMapper)
_titleActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
_titleActor.GetPositionCoordinate().SetValue(0.5, 0.95)

# set up the renderer and the render window
_ren = vtk.vtkRenderer()
_renWin = vtk.vtkRenderWindow()
_renWin.SetSize(640, 480)
_renWin.AddRenderer(_ren)

# add the actor
_ren.AddActor(_actor)
_ren.AddActor(_titleActor)

# render the scene
_iren = vtk.vtkRenderWindowInteractor()
_iren.SetRenderWindow(_renWin)
_iren.Initialize()
_renWin.Render()
_iren.Start()

# convert the render window to an image
_win2imgFilter = vtk.vtkWindowToImageFilter()
_win2imgFilter.SetInput(_renWin)

# save the image to file
_outWriter = vtk.vtkPNGWriter()
_outWriter.SetInput(_win2imgFilter.GetOutput())
_outWriter.SetFileName("ballPlot.png")
_outWriter.Write()

# vim: expandtab shiftwidth=4:

