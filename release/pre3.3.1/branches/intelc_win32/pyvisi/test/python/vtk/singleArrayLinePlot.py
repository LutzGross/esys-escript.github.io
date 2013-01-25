# $Id: singleArrayLinePlot.py,v 1.1 2005/05/05 01:57:53 paultcochrane Exp $

"""
Example of plotting a curve using only one input array with pyvisi 
"""

# set up some data to plot
from Numeric import *

x = arange(0,2*pi,0.1, typecode=Float)
y = sin(x)

#### original vtk code

import vtk

# set up the renderer and the render window
_ren = vtk.vtkRenderer()
_renWin = vtk.vtkRenderWindow()
_renWin.AddRenderer(_ren)

# do a quick check to make sure x and y are same length
if len(x) != len(y):
    raise ValueError, "x and y vectors must be same length"

# set up the x and y data arrays to be able to accept the data (code
# here adapted from the C++ of a vtk-users mailing list reply by Sander
# Niemeijer)
_xData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_xData.SetNumberOfTuples(len(y))

_yData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
_yData.SetNumberOfTuples(len(y))

# put the data into the data arrays
for i in range(len(x)):
    _xData.SetTuple1(i,i+1)
    _yData.SetTuple1(i,y[i])

# create a field data object 
# (I think this is as a containter to hold the data arrays)
_fieldData = vtk.vtkFieldData()
_fieldData.AllocateArrays(2)
_fieldData.AddArray(_xData)
_fieldData.AddArray(_yData)

# now put the field data object into a data object so that can add it as
# input to the xyPlotActor
_dataObject = vtk.vtkDataObject()
_dataObject.SetFieldData(_fieldData)

# set up the actor
_plot = vtk.vtkXYPlotActor()
_plot.AddDataObjectInput(_dataObject)

# set the title and stuff
_plot.SetTitle("Example 2D plot")
_plot.SetXTitle("index")
_plot.SetYTitle("y")
_plot.SetXValuesToValue()

# set which parts of the data object are to be used for which axis
_plot.SetDataObjectXComponent(0,0)
_plot.SetDataObjectYComponent(0,1)

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
_outWriter.SetFileName("singleArrayLinePlot.png")
_outWriter.Write()

# pause for input
raw_input('Press enter to continue...\n')


