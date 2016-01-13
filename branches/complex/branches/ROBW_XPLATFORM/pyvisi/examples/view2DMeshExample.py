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

# $Id: view2DMeshExample.py,v 1.3 2005/02/08 08:27:26 paultcochrane Exp $

## @file view2DMeshExample.py

"""
Example of loading a two-dimensional mesh into pyvisi and viewing it

Will hopefully help me write a decent interface.
""" 

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


