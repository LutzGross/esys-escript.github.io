# $Id: basicContourExample.py,v 1.7 2005/04/29 00:32:20 paultcochrane Exp $

"""
Example of a basic contour plot

Will hopefully help me write a decent interface.
"""

import sys,os

# import the python visualisation interface
from pyvisi import *

# original vtk code
import vtk

import numarray

# generate the x and y grid data
x = numarray.arrayrange(-1, 1, stride=0.5, type='Float')
y = numarray.arrayrange(-1, 1, stride=0.5, type='Float')

# generate a matrix of repeated x values (c.f. repmat() in matlab)
xm = numarray.zeros([len(x), len(y)], type='Float')
for i in range(len(y)):
    xm[:,i] = x

# generate a matrix of repeated y values (c.f. repmat() in matlab)
ym = numarray.zeros([len(x), len(y)], type='Float')
for i in range(len(x)):
    ym[i,:] = y

sigma = 0.2  # the spread of the distribution

# generate the distribution
distn = numarray.exp(-(xm*xm + ym*ym)/sigma)

# convert the x data into vtkFloatArray objects
xvtk = vtk.vtkFloatArray()
for i in range(len(x)):
    xvtk.InsertNextValue(x[i])

# convert the y data into vtkFloatArray objects
yvtk = vtk.vtkFloatArray()
for i in range(len(y)):
    yvtk.InsertNextValue(y[i])

# convert the distribution data into vtkFloatArray objects
distnvtk = vtk.vtkFloatArray()
for i in range(len(x)):
    for j in range(len(y)):
        distnvtk.InsertNextValue(distn[i,j])

# make the points to put into the grid
points = vtk.vtkPoints()
values = vtk.vtkFloatArray()
count = 0
for i in xrange(len(x)):
    for j in xrange(len(y)):
        points.InsertPoint(count, x[i], y[j], 0)
        values.InsertValue(count, distn[i,j])
        count += 1

# now make the strips (whatever they are...)
#strips = vtk.vtkCellArray()
#strips.InsertNextCell(len(x)*len(y))  # number of points
#for i in xrange(len(x)*len(y)):
    #strips.InsertCellPoint(i)

#strips.InsertCellPoint(0)
#strips.InsertCellPoint(1)
#strips.InsertCellPoint(7)
#strips.InsertCellPoint(6)
#strips.InsertCellPoint(2)
#strips.InsertCellPoint(3)
#strips.InsertCellPoint(5)
#strips.InsertCellPoint(4)

strips = vtk.vtkCellArray()
p2c = vtk.vtkPointDataToCellData()
p2c.SetInput(points)

# set up the polygonal data object
polyData = vtk.vtkPolyData()
polyData.SetPoints(points)
polyData.SetStrips(strips)
polyData.GetPointData().SetScalars(values)

warp = vtk.vtkWarpScalar()
warp.SetInput(polyData)
warp.SetScaleFactor(0.5)

contMapper = vtk.vtkPolyDataMapper()
contMapper.SetInput(warp.GetPolyDataOutput())
contMapper.SetScalarRange(polyData.GetScalarRange())

contActor = vtk.vtkActor()
contActor.SetMapper(contMapper)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.AddActor(contActor)
renWin.SetSize(400,400)
ren.SetBackground(1,1,1)
iren.Initialize()
renWin.Render()
iren.Start()
#raw_input("Press enter to continue")


# vim: expandtab shiftwidth=4:
