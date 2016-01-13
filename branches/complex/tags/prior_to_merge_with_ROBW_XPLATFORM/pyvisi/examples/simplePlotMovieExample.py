# $Id: simplePlotMovieExample.py,v 1.2 2005/03/07 04:17:40 paultcochrane Exp $

"""
Example of plotting a changing function with pyvisi 

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

x = arange(10, typecode=Float)
y = x**2

# plot it using one of the three methods
if method == 'pyvisi':

    # example code for how a user would write a script in pyvisi
    from esys.pyvisi import *          # base level visualisation stuff
    #from esys.pyvisi.utils import *   # pyvisi specific utils
    # import the objects to render the scene using the specific renderer
    from esys.pyvisi.renderers.gnuplot import *   # gnuplot
    #from esys.pyvisi.renderers.vtk import *       # vtk
    
    # define the scene object
    # a Scene is a container for all of the kinds of things you want to put 
    # into your plot for instance, images, meshes, arrow/vector/quiver plots, 
    # contour plots, spheres etc.
    scene = Scene()
    
    # create a LinePlot object
    plot = LinePlot(scene)
    
    # add some helpful info to the plot
    plot.title = 'Example 2D plot'
    plot.xlabel = 'x'
    plot.ylabel = 'x^2'

    plot.linestyle = 'lines'
    
    for i in range(100):
        # assign some data to the plot
        plot.setData(x, y)
    
        # render the scene to screen
        #scene.render(pause=True, interactive=True)
    
        # save the scene out to file
        scene.save(fname="simplePlotMovieExample%03d.png"%i, format=PngImage())
        #scene.save(fname="simplePlotMovieExample.ps", format=PsImage())
        y = y*0.9

elif method == 'gnuplot':
    #### original gnuplot code
    
    import Gnuplot

    # set the plot up
    _gnuplot = Gnuplot.Gnuplot()
    _gnuplot.title('Example 2D plot')
    _gnuplot.xlabel('x')
    _gnuplot.ylabel('x^2')

    for i in range(100):
        # set up the data
        _data = Gnuplot.Data(x, y, with='lines')
    
        # set up to save to file
        _gnuplot('set terminal png')
        _gnuplot('set output \"simplePlotMovieExample%03d.png\"'%i)
    
        # save it
        _gnuplot.plot(_data)

        y = y*0.9

    raw_input('Press enter to continue...\n')

elif method == 'vtk':
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
    _xData.SetNumberOfTuples(len(x))

    _yData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)
    _yData.SetNumberOfTuples(len(y))

    # put the data into the data arrays
    for i in range(len(x)):
        _xData.SetTuple1(i,x[i])
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
    _plot.SetXTitle("x")
    _plot.SetYTitle("x^2")
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
    _outWriter.SetFileName("simplePlotExample.png")
    _outWriter.Write()

    # pause for input
    #raw_input('Press enter to continue...\n')

else:
    print "Eeek!  What plotting method am I supposed to use???"

# vim: expandtab shiftwidth=4:

