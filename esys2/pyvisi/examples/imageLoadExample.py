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

# $Id$

## @file imageLoadExample.py

"""
Example of loading and viewing an image using pyvisi

Will hopefully help me write a decent interface.
"""

import sys,os

# this means that one can run the script from the examples directory
sys.path.append('../')

# import the python visualisation interface
from pyvisi import *

# set up a scene, using vtk to render it
scene = Scene(renderer="vtk")

# add a jpeg image to the scene, and then load the file
jpegImage = JpegImage(scene)
jpegImage.load(file="Flinders_eval.jpg")
jpegImage.render()  # this should be done at the scene.render step

# render the scene, pausing so that the opengl window doesn't disappear
scene.render(pause=True,interactive=True)

sys.exit()


# this is the original vtk code

import vtk 
ren = vtkRenderer()
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)

jpegReader = vtkJPEGReader()
jpegReader.SetFileName("Flinders_eval.jpg")

imgActor = vtkImageActor()
imgActor.SetInput(jpegReader.GetOutput())

ren.AddActor(imgActor)
renWin.SetSize(400,400)
ren.SetBackground(0.1,0.2,0.4)
renWin.Render()
raw_input("Press any key to continue")

# vim: expandtab shiftwidth=4:
