"""
View the Flinders ranges, with the fault lines lying on top

Will hopefully help me write a decent interface.

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

import sys,os

vtk = True

if not vtk:

    # this means that one can run the script from the examples directory
    sys.path.append('../')
    
    # import the python visualisation interface
    from esys.pyvisi import *
    
    # start a scene, using vtk as the renderer
    scene = Scene(renderer='vtk')
    scene.setBackgroundColor(1,1,1)

    # load some vtk data
    data = Data(scene)
    data.load(file="Flinders_ranges.vtk")

    # set up the plot
    plot = MeshPlot(scene)
    plot.setData(data)

    # add the plot to the scene
    scene.add(plot)
    
    # render the scene
    scene.render(pause=True,interactive=True)
    
    # the interactive flag means whether or not to set up and use
    # a window interactor object (if available)
	
    # put an exit in here so that we don't run the vtk code
    sys.exit()

else:

    # here is the original vtk code
    
    import vtk
    
    # set up the renderer and render window
    _renderer = vtk.vtkRenderer()
    _renderWindow = vtk.vtkRenderWindow()
    _renderWindow.AddRenderer(_renderer)
    _renderWindow.SetSize(640,480)
    _renderer.SetBackground(1,1,1)
    
    # load the vtk file
    _dataReader = vtk.vtkDataSetReader()
    _dataReader.SetFileName("Flinders_ranges.vtk")

    # load the jpeg image
    _jpegReader = vtk.vtkJPEGReader()
    _jpegReader.SetFileName("Flinders_eval.jpg")

    # map the image onto a plane
    _tex = vtk.vtkTexture()
    _tex.SetInput(_jpegReader.GetOutput())

    _plane = vtk.vtkPlaneSource()
    _planeMapper = vtk.vtkPolyDataMapper()
    _planeMapper.SetInput(_plane.GetOutput())
    _planeActor = vtk.vtkActor()
    _planeActor.SetMapper(_planeMapper)
    _planeActor.SetTexture(_tex)

    # alternatively, just grab the image
    _jpgImgActor = vtk.vtkImageActor()
    _jpgImgActor.SetInput(_jpegReader.GetOutput())

    # add the actor to the scene
    #_renderer.AddActor(_planeActor)
    #_renderer.AddActor(_jpgImgActor)

    # set up the data
    _dataMapper = vtk.vtkDataSetMapper()
    _dataMapper.SetInput(_dataReader.GetOutput())
    #_dataMapper.ScalarVisibilityOff()  # what exactly does this do?

    # set up the actor for viewing it
    _dataActor = vtk.vtkActor()
    _dataActor.SetMapper(_dataMapper)
    _dataActor.GetProperty().SetColor(0.2,0.2,0.2)
    _dataActor.GetProperty().SetRepresentationToWireframe()
    
    # add the actor to the scene
    #_renderer.AddActor(_dataActor)

    # try with an assembly
    _assembly = vtk.vtkAssembly()
    _assembly.AddPart(_planeActor)
    _assembly.AddPart(_dataActor)

    # add the assembly to the scene
    _renderer.AddActor(_assembly)

    # now see what was produced, with interactive playing stuff
    _iRenderer = vtk.vtkRenderWindowInteractor()
    _iRenderer.SetRenderWindow(_renderWindow)
    _iRenderer.Initialize()
    _renderWindow.Render()
    _iRenderer.Start()

    raw_input("Press enter to continue")

# vim: expandtab shiftwidth=4:
