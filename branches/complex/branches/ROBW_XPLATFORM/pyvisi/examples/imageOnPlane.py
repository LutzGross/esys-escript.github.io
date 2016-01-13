# Copyright (C) 2004 Paul Cochrane
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# $Id: imageOnPlane.py,v 1.2 2005/01/11 05:47:51 paultcochrane Exp $

## @file imageLoadExample.py

"""
Example of loading an image into pyvisi and mapping onto a plane

Will hopefully help me write a decent interface.
""" 

import sys,os

# this means that one can run the script from the examples directory
sys.path.append('../')

# import the python visualisation interface
from pyvisi import *
# import the vtk stuff
from pyvisi.renderers.vtk import *

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


