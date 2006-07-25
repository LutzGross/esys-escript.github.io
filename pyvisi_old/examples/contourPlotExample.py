# $Id: contourPlotExample.py,v 1.8 2005/05/24 01:32:12 paultcochrane Exp $

"""
Example of contour plotting with pyvisi 

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


# what plotting method are we using?
method = 'pyvisi'

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

# plot it with either gnuplot, vtk or pyvisi
if method == 'pyvisi':
    #### pyvisi version of code

    # import the general pyvisi stuff
    from esys.pyvisi import *
    # import the gnuplot overrides of the interface
    from esys.pyvisi.renderers.gnuplot import *
    #from esys.pyvisi.renderers.vtk import *

    # define a scene object
    # a Scene is a container for all of the kinds of things you want to put
    # into your plot, for instance, images, meshes, arrow/vector/quiver
    # plots, contour plots, spheres etc.
    scene = Scene()

    # create a ContourPlot object
    plot = ContourPlot(scene)

    # add some helpful info to the plot
    plot.title = 'Example contour plot'
    plot.xlabel = 'x'
    plot.ylabel = 'y'

    # assign the data to the plot
    # this version assumes that we have x, then y, then z and that z is 2D
    # and that x and y are 1D arrays
    plot.setData(x,y,z)
    # alternative syntax
    #plot.setData(xData=x, yData=y, zData=z)
    # or (but more confusing depending upon one's naming conventions)
    #plot.setData(x=x, y=y, z=z)

    # render the scene to screen
    scene.render(pause=True, interactive=True)

    # save the scene to file
    plot.setData(x,y,z)  # have to do this now because we've already
                         # render()ed the scene.  This requirement will be
                         # removed in the future
    scene.save(fname="contourPlotExample.png", format=PngImage())

elif method == 'vtk':
    #### original vtk code

    import vtk

    # set up the data
    _points = vtk.vtkPoints()
    _scalars = vtk.vtkFloatArray()

    _index = 0
    for _i in range(len(x)):
        for _j in range(len(y)):
            _points.InsertPoint(_index, x[_i], y[_j], 0)
            _scalars.InsertValue(_index, z[_i,_j])

    zMin = min(min(z))
    zMax = max(max(z))

    _data = vtk.vtkUnstructuredGrid()
    _data.SetPoints(_points)
    _data.GetPointData().SetScalars(_scalars);

    # set up the contour
    _plotContour = vtk.vtkContourGrid()
    _plotContour.SetInput(_data)
    _plotContour.GenerateValues(5, zMin, zMax)

    # set up the mapper
    _plotMapper = vtk.vtkPolyDataMapper()
    _plotMapper.SetInput(_plotContour.GetOutput())
    _plotMapper.SetScalarRange(zMin, zMax)

    # set up the actor
    _plotActor = vtk.vtkActor()
    _plotActor.SetMapper(_plotMapper)

    # use a scalar bar
    #_scalarBar = vtk.vtkScalarBarActor()

    #_lut = vtk.vtkLookupTable()
    #_planeMapper = vtk.vtkPolyDataMapper()
    #_planeMapper.SetInput(_data)
    #_planeMapper.SetLookupTable(_lut)
    #_planeMapper.SetScalarRange(0.0, 1.0)
    #_scalarBarActor = vtk.vtkActor()
    #_scalarBarActor.SetMapper(_planeMapper)

    # set up the renderer and the render window
    _ren = vtk.vtkRenderer()
    _ren.SetBackground(1, 1, 1)
    _renWin = vtk.vtkRenderWindow()
    _renWin.AddRenderer(_ren)

    # add the actor to the renderer
    _ren.AddActor(_plotActor)
    #_ren.AddActor(_scalarBarActor)

    # render the scene
    _iren = vtk.vtkRenderWindowInteractor()
    _iren.SetRenderWindow(_renWin)
    _iren.Initialize()
    _renWin.Render()
    _iren.Start()

    #raw_input('Press enter to continue...')

elif method == 'plplot':
    pass
else:
    print "Eeek!  What plotting method am I supposed to use???"

# vim: expandtab shiftwidth=4:
