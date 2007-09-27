"""
Example of loading a two-dimensional mesh into pyvisi and viewing it

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

# this means that one can run the script from the examples directory
sys.path.append('../')

# import the python visualisation interface
from esys.pyvisi import *
# import vtk stuff
from esys.pyvisi.renderers.vtk import *

# start a scene, using vtk as the renderer
scene = Scene()
scene.setSize()
scene.setBackgroundColor(0,0,0)

# render the scene
scene.render(pause=True)

# put an exit in here so that we don't run the vtk code
sys.exit()

# here is the original vtk code

import vtk

# set up the renderer and render window
_renderer = vtk.vtkRenderer()
_renderWindow = vtk.vtkRenderWindow()
_renderWindow.AddRenderer(_renderer)
_renderWindow.SetSize(1024,768)
_renderer.SetBackground(0,0,0)

# load the jpeg file
_jpegReader = vtk.vtkJPEGReader()
_jpegReader.SetFileName("Flinders_eval.jpg")

# set the jpeg file to the texture
_tex = vtk.vtkTexture()
_tex.SetInput(_jpegReader.GetOutput())

# get the plane and add it to the "scene"
_plane = vtk.vtkPlaneSource()
_planeMapper = vtk.vtkPolyDataMapper()
_planeMapper.SetInput(_plane.GetOutput())
_planeActor = vtk.vtkActor()
_planeActor.SetMapper(_planeMapper)
_planeActor.SetTexture(_tex)
_renderer.AddActor(_planeActor)

# set up a camera
_camera = vtk.vtkCamera()
_camera.SetFocalPoint(0,0,0)
_camera.SetPosition(1,1,1)
_camera.ComputeViewPlaneNormal()
_camera.SetViewUp(1,0,0)
_camera.OrthogonalizeViewUp()

_renderer.SetActiveCamera(_camera)
_renderer.ResetCamera()

# now see what was produced
_renderWindow.Render()

raw_input("Press enter to continue")


