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

## @file plane.py

"""
The classes associated with Planes
"""

from common import _debug
from item import Item

class Plane(Item):
    """
    Generic class for Plane objects
    """

    def __init__(self,scene):
        """
        Initialisation of the Plane object
        """
        if _debug: print "\tCalled Plane.__init__()"

        self.renderer = scene.renderer
        return

    def mapImageToPlane(self,image):
        """
        Maps an Image object onto a Plane object
        """
        if _debug: print "\tCalled Plane.mapImageToPlane()"

        # need to work out the name of the internal image object name
        imgObjectName = image.readerName

        self.renderer.addToEvalStack("# Plane.mapImageToPlane()\n")
        self.renderer.addToEvalStack("_tex = vtk.vtkTexture()\n")
        evalString = "_tex.SetInput(%s.GetOutput())\n" % imgObjectName
        self.renderer.addToEvalStack(evalString)
        self.renderer.addToEvalStack("_plane = vtk.vtkPlaneSource()\n")
        self.renderer.addToEvalStack(\
                "_planeMapper = vtk.vtkPolyDataMapper()\n")
        self.renderer.addToEvalStack(\
                "_planeMapper.SetInput(_plane.GetOutput())\n")
        return

    def render(self):
        """
        Perform Plane object specific (pre)rendering tasks
        """
        if _debug: print "\tCalled Plane.mapImageToPlane()"

        self.renderer.addToEvalStack("# Plane.render()\n")
        self.renderer.addToEvalStack("_planeActor = vtk.vtkActor()\n")
        self.renderer.addToEvalStack("_planeActor.SetMapper(_planeMapper)\n")
        self.renderer.addToEvalStack("_planeActor.SetTexture(_tex)\n")
        self.renderer.addToEvalStack("_renderer.AddActor(_planeActor)\n")

        return
    
# vim: expandtab shiftwidth=4:
