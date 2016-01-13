# $Id: offsetLinePlot.py,v 1.1 2005/11/08 06:51:39 paultcochrane Exp $ 
"""
Example of plotting multiple curves offset from each other

This is especially handy for people plotting seismic data
"""

# set up some data to plot
from Numeric import *

x = arange(0,2*pi,0.01, typecode=Float)
y1 = sin(x)
y2 = cos(x)
y3 = cos(x)**2
y4 = sin(2*x)
y5 = cos(3*x)
y6 = sin(20*x)

#### original vtk code

import vtk

# set up the renderer and the render window
_ren = vtk.vtkRenderer()
_renWin = vtk.vtkRenderWindow()
_renWin.AddRenderer(_ren)

# do a quick check to make sure x and y are same length
if len(x) != len(y1):
    raise ValueError, "x and y vectors must be same length"

if len(x) != len(y2):
    raise ValueError, "x and y vectors must be same length"

if len(x) != len(y3):
    raise ValueError, "x and y vectors must be same length"

if len(x) != len(y4):
    raise ValueError, "x and y vectors must be same length"

if len(x) != len(y5):
    raise ValueError, "x and y vectors must be same length"

if len(x) != len(y6):
    raise ValueError, "x and y vectors must be same length"

### set up the data to be offset
# concatenate the data
yAll = concatenate( [y1, y2, y3, y4, y5, y6] )

yMax = max(yAll)
yMin = min(yAll)

# keep the data apart a bit
const = 0.1*(yMax-yMin)

# don't need to worry about y1: it's the first data series
shift = yMax - yMin + const
y2 = y2 + shift
y3 = y3 + 2*shift
y4 = y4 + 3*shift
y5 = y5 + 4*shift
y6 = y6 + 5*shift

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

_yData4 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_yData4.SetNumberOfTuples(len(y4))

_yData5 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_yData5.SetNumberOfTuples(len(y5))

_yData6 = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_yData6.SetNumberOfTuples(len(y6))

# put the data into the data arrays
for i in range(len(x)):
    _xData.SetTuple1(i,x[i])
    _yData1.SetTuple1(i,y1[i])
    _yData2.SetTuple1(i,y2[i])
    _yData3.SetTuple1(i,y3[i])
    _yData4.SetTuple1(i,y4[i])
    _yData5.SetTuple1(i,y5[i])
    _yData6.SetTuple1(i,y6[i])

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

_fieldData4 = vtk.vtkFieldData()
_fieldData4.AllocateArrays(2)
_fieldData4.AddArray(_xData)
_fieldData4.AddArray(_yData4)

_fieldData5 = vtk.vtkFieldData()
_fieldData5.AllocateArrays(2)
_fieldData5.AddArray(_xData)
_fieldData5.AddArray(_yData5)

_fieldData6 = vtk.vtkFieldData()
_fieldData6.AllocateArrays(2)
_fieldData6.AddArray(_xData)
_fieldData6.AddArray(_yData6)

# now put the field data object into a data object so that can add it as
# input to the xyPlotActor
_dataObject1 = vtk.vtkDataObject()
_dataObject1.SetFieldData(_fieldData1)

_dataObject2 = vtk.vtkDataObject()
_dataObject2.SetFieldData(_fieldData2)

_dataObject3 = vtk.vtkDataObject()
_dataObject3.SetFieldData(_fieldData3)

_dataObject4 = vtk.vtkDataObject()
_dataObject4.SetFieldData(_fieldData4)

_dataObject5 = vtk.vtkDataObject()
_dataObject5.SetFieldData(_fieldData5)

_dataObject6 = vtk.vtkDataObject()
_dataObject6.SetFieldData(_fieldData6)

# set up the actor
_plot = vtk.vtkXYPlotActor()
_plot.AddDataObjectInput(_dataObject1)
_plot.AddDataObjectInput(_dataObject2)
_plot.AddDataObjectInput(_dataObject3)
_plot.AddDataObjectInput(_dataObject4)
_plot.AddDataObjectInput(_dataObject5)
_plot.AddDataObjectInput(_dataObject6)

# set the title and stuff
_plot.SetTitle("Example 2D plot with offsets")
_plot.SetXTitle("x")
_plot.SetYTitle("y")
_plot.SetXValuesToValue()

# set which parts of the data object are to be used for which axis
_plot.SetDataObjectXComponent(0,0)
_plot.SetDataObjectYComponent(0,1)
_plot.SetDataObjectYComponent(1,1)
_plot.SetDataObjectYComponent(2,1)
_plot.SetDataObjectYComponent(3,1)
_plot.SetDataObjectYComponent(4,1)
_plot.SetDataObjectYComponent(5,1)

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
_outWriter.SetFileName("offsetLinePlot.png")
_outWriter.Write()

# pause for input
#raw_input('Press enter to continue...\n')

# vim: expandtab shiftwidth=4:
