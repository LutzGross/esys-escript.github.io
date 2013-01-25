"""
Example of loading an image into pyvisi and mapping onto a plane

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
from pyvisi_old import *
# import the vtk stuff
from pyvisi_old.renderers.vtk import *

# start a scene, using vtk as the renderer
scene = Scene()

# load a jpeg image into the scene
jpegImage = JpegImage(scene)
jpegImage.load(file="Flinders_eval.jpg")

# define a plane to map image to, and map it
plane = Plane(scene)
plane.mapImageToPlane(jpegImage)
plane.render()

# maybe instead of this render() mechanism, I should be explictly adding
# objects to be rendered at the end.  Also, instead of trying to (say) render
# absolutely everything in that is intantiated, maybe have the add() mechanism
# that Alexei recommended.  Say instead, jpegImage.addToScene(scene), which
# would have the ability to add the same object to different scenes....  And
# then all "added" objects would be rendered in the end.  Need to think on
# this one...

# initialise a camera
camera = Camera(scene)

# I would far rather use this code, but atm, I can't work out how to do it
# scene.render()
# so I'm going to use this
#for i in range(360):
    #scene.vtkCommand("_renderWindow.Render()")
    #camera.setElevation(1)
    #camera.setAzimuth(1)
    #scene.vtkCommand("_renderer.ResetCamera()")
    #evalStr = "print \"Angle is %d\"" % i
    #scene.vtkCommand(evalStr)

scene.setBackgroundColor(0,0,0)

# render the scene
scene.render(pause=True,interactive=True)

# put an exit in here so that we don't run the vtk code
sys.exit()

# here is the original vtk code

import vtk

# set up the renderer and render window
_renderer = vtk.vtkRenderer()
_renderWindow = vtk.vtkRenderWindow()
_renderWindow.AddRenderer(_renderer)
_renderWindow.SetSize(1024,768)
_renderer.SetBackground(1,1,1)

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


