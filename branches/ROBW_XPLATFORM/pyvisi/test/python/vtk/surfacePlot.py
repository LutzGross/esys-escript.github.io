# $Id: surfacePlot.py,v 1.1 2005/11/08 07:26:31 paultcochrane Exp $

"""
Example of plotting surfaces
"""

# set up some data to plot
from Numeric import *

# the x and y axes
x = arange(-2,2,0.2, typecode=Float)
y = arange(-2,3,0.2, typecode=Float)

# pick some interesting function to generate the data in the third dimension
# this is the one used in the matlab docs: z = x*exp(-x^2-y^2)
z = zeros((len(x),len(y)), typecode=Float)

# boy do *I* feel old fashioned writing it this way
# surely there's another way to do it: - something to do later
for i in range(len(x)):
    for j in range(len(y)):
	z[i,j] = x[i]*exp(-x[i]*x[i] - y[j]*y[j])

#### original vtk code

import vtk

# set up the points
_points = vtk.vtkPoints()
_points.SetNumberOfPoints(len(x))
_index = 0
for i in range(len(x)):
    for j in range(len(y)):
	_points.InsertPoint(_index, x[i], y[j], 0)
	_index += 1

# create the data
_data = vtk.vtkFloatArray()
_data.SetNumberOfComponents(1)
_data.SetNumberOfValues(len(x)*len(y))
_index = 0
for i in range(len(x)):
    for j in range(len(y)):
	_data.InsertValue(_index, z[i][j])
	_index += 1

# set up the grid (it's polydata since we're doing a Delaunay2D)
_grid = vtk.vtkPolyData()
_grid.SetPoints(_points)
_grid.GetPointData().SetScalars(_data)

# calculate the min and max of the z data
_zMin = min(z.flat)
_zMax = max(z.flat)

# make a lookup table for the colour map and invert it (colours look
# better when it's inverted)
_lut = vtk.vtkLookupTable()
_refLut = vtk.vtkLookupTable()
_lut.Build()
_refLut.Build()
for i in range(256):
    _lut.SetTableValue(i, _refLut.GetTableValue(255-i))

# triangulate them
_delaunay = vtk.vtkDelaunay2D()
_delaunay.SetInput(_grid)
_delaunay.SetTolerance(0.001)

# warp the surface to generate the surface or "carpet"
_warp = vtk.vtkWarpScalar()
_warp.SetInput(_delaunay.GetOutput())
_warp.SetScaleFactor(1.0)

# set up the mapper
_mapper = vtk.vtkPolyDataMapper()
_mapper.SetInput(_warp.GetOutput())
_mapper.SetLookupTable(_lut)
_mapper.SetScalarRange(_zMin, _zMax)

# set up the actor
_plotActor = vtk.vtkActor()
_plotActor.SetMapper(_mapper)

# set up the text properties for nice text
_font_size = 20
_textProp = vtk.vtkTextProperty()
_textProp.SetFontSize(_font_size)
_textProp.SetFontFamilyToArial()
_textProp.BoldOff()
_textProp.ItalicOff()
_textProp.ShadowOff()

# use a scalar bar
_scalarBar = vtk.vtkScalarBarActor()
_scalarBar.SetLookupTable(_lut)
_scalarBar.SetWidth(0.1)
_scalarBar.SetHeight(0.7)
_scalarBar.SetPosition(0.9, 0.2)
_scalarBar.SetTitle("z")  # don't need this I don't think

# make a title
_title = vtk.vtkTextMapper()
_title.SetInput("Example surface plot")

# make the title text use the text properties
_titleProp = _title.GetTextProperty()
_titleProp.ShallowCopy(_textProp)
_titleProp.SetJustificationToCentered()
_titleProp.SetVerticalJustificationToTop()

# make the actor for the title
_titleActor = vtk.vtkTextActor()
_titleActor.SetMapper(_title)
_titleActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
_titleActor.GetPositionCoordinate().SetValue(0.5, 0.95)

# set up the renderer and the render window
_ren = vtk.vtkRenderer()
_renWin = vtk.vtkRenderWindow()
_renWin.AddRenderer(_ren)
_renWin.SetSize(640,480)

# add the actor to the renderer
_ren.AddActor(_plotActor)
_ren.AddActor(_scalarBar)
_ren.AddActor(_titleActor)

# set up the camera, and associate the text with the camera
_camera = _ren.GetActiveCamera()
_gridCentre = _grid.GetCenter()
_camera.SetFocalPoint(_gridCentre)
_camPos = _camera.GetPosition()
_pos = (_camPos[0] - _gridCentre[0], 
	_camPos[1] - _gridCentre[1],
	_camPos[2] - _gridCentre[2])
_radius = sqrt(_pos[0]*_pos[0] + _pos[1]*_pos[1] + _pos[2]*_pos[2])

# these are the matlab defaults...  they look ok, so use them
# note that azimuth and elevation mean slightly different things in vtk
# to what they mean in matlab: they aren't absolute, they are relative
# to the current position
_azimuth = -37.5
_elevation = 30

# take the position where I am as being 0 azimuth, and 90 elevation

# position the camera appropriately
_xPos = _gridCentre[0] + _radius*sin(math.pi*_azimuth/180.0)
_yPos = _gridCentre[1] - _radius*cos(math.pi*_elevation/180.0) \
	+ _radius*cos(math.pi*_azimuth/180.0)
_zPos = _gridCentre[2] + _radius - _radius*sin(math.pi*_elevation/180.0)

_camera.SetPosition(_xPos, _yPos, _zPos) 

_camera.SetViewUp(0,0,1)

_ren.SetActiveCamera(_camera)

_ren.ResetCameraClippingRange()

# add some axes
_axes = vtk.vtkCubeAxesActor2D()
_axes.SetInput(_warp.GetOutput())
_axes.SetCamera(_ren.GetActiveCamera())
_axes.SetLabelFormat("%6.2g")
_axes.SetFlyModeToOuterEdges()
_axes.SetFontFactor(0.8)
_axes.SetAxisTitleTextProperty(_textProp)
_axes.SetAxisLabelTextProperty(_textProp)
_axes.SetXLabel("x")
_axes.SetYLabel("y")
_axes.SetZLabel("z")
_axes.SetNumberOfLabels(5)
_ren.AddProp(_axes)

# render the scene
_iren = vtk.vtkRenderWindowInteractor()
_iren.SetRenderWindow(_renWin)
_iren.Initialize()
_renWin.Render()
_iren.Start()

# vim: expandtab shiftwidth=4:

