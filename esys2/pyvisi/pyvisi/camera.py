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

## @file camera.py

"""
@brief Class and functions associated with a pyvisi Camera object
"""

from common import _debug
from item import Item

class Camera(Item):
    """
    @brief Camera class
    """
    def __init__(self,scene):
        """
        Initialisation of the Camera object

        @param scene The Scene object to add the Camera object to
        """
        if _debug: print "\tCalled Camera.__init__()"

        # default x,y,z positions of Camera (specific to vtk)
        self.xPos = 0.0
        self.yPos = 0.0
        self.zPos = 3.0

        # default x,y,z positions of the Camers's focal point (specific to vtk)
        self.xFocalPoint = 0.0
        self.yFocalPoint = 0.0
        self.zFocalPoint = 0.0

        # keep a reference to the renderer so we can send stuff to it
        self.renderer = scene.renderer

        # some vtk initialisation commands
        self.renderer.addToEvalStack("# Camera.__init__()\n")
        self.renderer.addToEvalStack("_camera = _renderer.GetActiveCamera()\n")

        # initialise the position of the Camera
        self.setPosition(self.xPos, self.yPos, self.zPos)
        self.setFocalPoint(self.xFocalPoint, self.yFocalPoint, self.zFocalPoint)
        self.renderer.addToEvalStack("_renderer.SetActiveCamera(_camera)\n")
        # dunno if this next line is exactly necessary
        self.renderer.addToEvalStack("_renderer.ResetCamera()\n")
        return

    def setPosition(self,*pos):
        """
        Set position of camera within scene

        @param pos Position to set camera in terms of x,y,z coordinates
        """
        if _debug: print "\tCalled Camera.setPosition()"

        # I need to do some mucking around in here with coordinate systems
        # and so on, but at present, we'll just use vtk's coord system
        
        self.xPos = pos[0]
        self.yPos = pos[1]
        self.zPos = pos[2]

        # now to set the position
        self.renderer.addToEvalStack("# Camera.setPosition()\n")
        evalString = "_camera.SetPosition(%f,%f,%f)\n" % \
                (self.xPos, self.yPos, self.zPos)
        self.renderer.addToEvalStack(evalString)

        return

    def getPosition(self):
        """
        Get the position of Camera within Scene

        Returns the position in a tuple of form (xPos, yPos, zPos)
        """
        if _debug: print "\tCalled Camera.getPosition()"

        return (self.xPos, self.yPos, self.zPos)

    def setFocalPoint(self,*pos):
        """
        Sets the focal point of the Camera with the Scene

        @param pos Position to set the focal point
        """
        if _debug: print "\tCalled Camera.setFocalPoint()"

        # I need to do some mucking around in here with coordinate systems
        # and so on, but at present, we'll just use vtk's coord system
        
        self.xFocalPoint = pos[0]
        self.yFocalPoint = pos[1]
        self.zFocalPoint = pos[2]

        # now set the focal point position
        self.renderer.addToEvalStack("#Camera.setFocalPoint()\n")
        evalString = "_camera.SetFocalPoint(%f,%f,%f)\n" % \
                (self.xFocalPoint, self.yFocalPoint, self.zFocalPoint)
        self.renderer.addToEvalStack(evalString)

        return

    def getFocalPoint(self):
        """
        Get the position of the focal point of the Camera

        Returns the position of the focal point in a tuple of form 
        (xPos, yPos, zPos)
        """
        if _debug: print "\tCalled Camera.getFocalPoint()"

        return (self.xFocalPoint, self.yFocalPoint, self.zFocalPoint)

    def setElevation(self,elevation):
        """
        Set the elevation angle (in degrees) of the Camera

        @param elevation The elevation angle (in degrees) of the Camera
        """
        if _debug: print "\tCalled Camera.setElevation()"

        self.elevation = elevation
        evalString = "_camera.Elevation(%f)\n" % elevation
        self.renderer.addToEvalStack(evalString)

        return

    def getElevation(self):
        """
        Gets the elevation angle (in degrees) of the Camera
        """
        if _debug: print "\tCalled Camera.getElevation()"
        
        return self.elevation

    def setAzimuth(self,azimuth):
        """
        Set the azimuthal angle (in degrees) of the Camera

        @param azimuth The azimuthal angle (in degrees) of the Camera
        """
        if _debug: print "\tCalled Camera.setAzimuth()"

        self.azimuth = azimuth
        evalString = "_camera.Azimuth(%f)\n" % azimuth
        self.renderer.addToEvalStack(evalString)

        return

    def getAzimuth(self):
        """
        Get the azimuthal angle (in degrees) of the Camera
        """
        if _deubg: print "\tCalled Camera.getAzimuth()"

        return self.azimuth


# vim: expandtab shiftwidth=4:
