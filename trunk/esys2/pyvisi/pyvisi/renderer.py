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

## @file render.py

"""
This is the file for the Renderer class
"""

from common import _debug

class Renderer(object):
    """
    A generic object holding a renderer of a Scene().
    """

    def __init__(self,rendererName):
        """
        Initialisation of Renderer() class
    
        @param rendererName The name of the renderer backend to use
        """
        if _debug: print "\tCalled Renderer.__init__()"
        self.rendererName = rendererName
        if _debug: print "\tRenderer name is %s" % rendererName

        # initialise some attributes
        self.renderWindowWidth = 640
        self.renderWindowHeight = 480

        # initialise the evalstack
        self._evalStack = ""

        # initialise the renderer module
        self._initRendererModule(self.rendererName)

    def _initRendererModule(self,rendererName):
        """
        Initialise the renderer module

        @param rendererName The name of the renderer to initialise
        """
        if _debug: print "\tCalled Renderer._initRendererModule()"
        # do a check to see if the renderer used is one of the known renderers
        if rendererName == "vtk":
            self._evalStack += "# Renderer._initRendererModule\n"
            self._evalStack += "import vtk\n"
            self._evalStack += "_renderer = vtk.vtkRenderer()\n"
            # this next line should only be done if we have an active display
            self._evalStack += "_renderWindow = vtk.vtkRenderWindow()\n"
            self._evalStack += "_renderWindow.AddRenderer(_renderer)\n"
            self._evalStack += "_renderWindow.SetSize(%d,%d)\n" % \
                    (self.renderWindowWidth,self.renderWindowHeight)
        else:
            print "Unknown renderer"
            return None

        return

    def setRenderWindowWidth(self,width):
        """
        Sets the render window width
        
        @param width The width of the render window
        """
        self.renderWindowWidth = width
        return

    def setRenderWindowHeight(self,height):
        """
        Sets the render window height

        @param height The height of the render window
        """
        self.renderWindowHeight = height
        return

    def getRenderWindowWidth(self):
        """
        Gets the render window width
        """
        return self.renderWindowWidth

    def getRenderWindowHeight(self):
        """
        Gets the render window height
        """
        return self.renderWindowHeight

    def getEvalStack(self):
        """
        Gets the evaluation stack as it currently stands
        """
        return self._evalStack

    def addToEvalStack(self,evalString):
        """
        Method to add commands to the evaluation stack
        
        @param renderer The renderer whose evalStack to add to
        @param evalString The string of commands to be added to the evalStack
        """
        if _debug: print "\tCalled Renderer.addToEvalStack()"
        self._evalStack += evalString
        return

    def resetEvalStack(self):
        """
        Reset/flush the evaluation stack
        """
        if _debug: print "\tCalled Renderer.resetEvalStack()"
        self._evalStack = ""
        return

# vim: expandtab shiftwidth=4:
