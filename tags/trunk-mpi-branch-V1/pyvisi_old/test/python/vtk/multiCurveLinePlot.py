# $Id: multiCurveLinePlot.py,v 1.1 2005/05/05 01:57:53 paultcochrane Exp $

"""
Example of plotting multiple curves with pyvisi 
"""

# set up some data to plot
from Numeric import *

x = arange(0, 2*pi, 0.1, typecode=Float)
y1 = sin(x)
y2 = cos(x)
y3 = cos(x)**2

#### original vtk code

import vtk

# set up the renderer and the render window
_ren = vtk.vtkRenderer()
_renWin = vtk.vtkRenderWindow()
_renWin.AddRenderer(_ren)
_renWin.SetSize(640, 480)

# do a quick check to make sure x and y are same length
if len(x) != len(y1):
    raise ValueError, "x and y vectors must be same length"

if len(x) != len(y2):
    raise ValueError, "x and y vectors must be same length"

# set up the x and y data arrays to be able to accept the data (code
# here adapted from the C++ of a vtk-users mailing list reply by Sander
# Niemeijer)
_xData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_xData.SetNumberOfTuples(len(x))

_yData1 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_yData1.SetNumberOfTuples(len(y1))

_yData2 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_yData2.SetNumberOfTuples(len(y2))

_yData3 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_yData3.SetNumberOfTuples(len(y3))

# put the data into the data arrays
for i in range(len(x)):
    _xData.SetTuple1(i,x[i])
    _yData1.SetTuple1(i,y1[i])
    _yData2.SetTuple1(i,y2[i])
    _yData3.SetTuple1(i,y3[i])

# create a field data object 
# (I think this is as a containter to hold the data arrays)
_fieldData1 = vtk.vtkFieldData()
_fieldData1.AllocateArrays(2)
_fieldData1.AddArray(_xData)
_fieldData1.AddArray(_yData1)

_fieldData2 = vtk.vtkFieldData()
_fieldData2.AllocateArrays(2)
_fieldData2.AddArray(_xData)
_fieldData2.AddArray(_yData2)

_fieldData3 = vtk.vtkFieldData()
_fieldData3.AllocateArrays(2)
_fieldData3.AddArray(_xData)
_fieldData3.AddArray(_yData3)

# now put the field data object into a data object so that can add it as
# input to the xyPlotActor
_dataObject1 = vtk.vtkDataObject()
_dataObject1.SetFieldData(_fieldData1)

_dataObject2 = vtk.vtkDataObject()
_dataObject2.SetFieldData(_fieldData2)

_dataObject3 = vtk.vtkDataObject()
_dataObject3.SetFieldData(_fieldData3)

# set up the actor
_plot = vtk.vtkXYPlotActor()
_plot.AddDataObjectInput(_dataObject1)
_plot.AddDataObjectInput(_dataObject2)
_plot.AddDataObjectInput(_dataObject3)

# the size of the actor should be 80% of the render window
_plot.SetPosition(0.1, 0.1)  # (0.1 = (1.0 - 0.8)/2)
_plot.SetWidth(0.8)
_plot.SetHeight(0.8)

# set the title and stuff
_plot.SetTitle("Example 2D plot")
_plot.SetXTitle("x")
_plot.SetYTitle("y")
_plot.SetXValuesToValue()

 # set which parts of the data object are to be used for which axis
_plot.SetDataObjectXComponent(0, 0)
_plot.SetDataObjectYComponent(0, 1)
_plot.SetDataObjectYComponent(1, 1)
_plot.SetDataObjectYComponent(2, 1)

# set up the lookup table for the appropriate range of colours
_lut = vtk.vtkLookupTable()
_lut.Build()
_colours = []
for i in range(3):
    _colours.append(_lut.GetColor(float(i)/float(3-1)))

# change the colour of the separate lines
_plot.SetPlotColor(0, _colours[0][0], _colours[0][1], _colours[0][2])
_plot.SetPlotColor(1, _colours[1][0], _colours[1][1], _colours[1][2])
_plot.SetPlotColor(2, _colours[2][0], _colours[2][1], _colours[2][2])

_plot.GetXAxisActor2D().GetProperty().SetColor(0,0,0)
_plot.GetYAxisActor2D().GetProperty().SetColor(0,0,0)
_ren.SetBackground(1,1,1)

# add the actor
_ren.AddActor2D(_plot)

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
_outWriter.SetFileName("multiCurveLinePlot.png")
_outWriter.Write()

# pause for input
raw_input('Press enter to continue...\n')

