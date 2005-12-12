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

# $Id: vtkDataLoadExample.py,v 1.3 2005/11/08 05:13:04 paultcochrane Exp $

## @file view2DMeshExample.py

"""
Example of loading a two-dimensional mesh into pyvisi and viewing it

Will hopefully help me write a decent interface.
""" 

import sys,os

vtk = True

if not vtk:

    # this means that one can run the script from the examples directory
    sys.path.append('../')
    
    # import the python visualisation interface
    from pyvisi import *
    # import vtk stuff
    from pyvisi.renderers.vtk import *
    
    # start a scene, using vtk as the renderer
    scene = Scene()
    scene.setBackgroundColor(1,1,1)

    # load some vtk data
    data = Data(scene)
    data.load(file="Flinders_ranges.vtk")

    # set up the plot
    plot = MeshPlot(scene)
    plot.setData(data)

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
    #_dataReader.SetFileName("Flinders_ranges.vtk")
    _dataReader.SetFileName("t.vtk")
	
    # set up the data
    _dataMapper = vtk.vtkDataSetMapper()
    _dataMapper.SetInput(_dataReader.GetOutput())
    _dataMapper.ScalarVisibilityOff()  # what exactly does this do?

    # set up the actor for viewing it
    _dataActor = vtk.vtkActor()
    _dataActor.SetMapper(_dataMapper)
    _dataActor.GetProperty().SetColor(0.2,0.2,0.2)
    _dataActor.GetProperty().SetRepresentationToWireframe()
    
    # add the actor to the scene
    _renderer.AddActor(_dataActor)

    # now see what was produced, with interactive playing stuff
    _iRenderer = vtk.vtkRenderWindowInteractor()
    _iRenderer.SetRenderWindow(_renderWindow)
    _iRenderer.Initialize()
    _renderWindow.Render()
    _iRenderer.Start()

    raw_input("Press enter to continue")

# vim: expandtab shiftwidth=4:
