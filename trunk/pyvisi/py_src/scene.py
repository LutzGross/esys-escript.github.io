# Copyright (C) 2004-2005 Paul Cochrane
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

# $Id: scene.py,v 1.23 2005/11/03 05:55:59 paultcochrane Exp $

## @file scene.py

"""
Class and functions associated with a pyvisi Scene
"""

# generic imports
from pyvisi.common import debugMsg, overrideWarning, fileCheck

from pyvisi.renderer import Renderer

__revision__ = '$Revision: 1.23 $'

class Scene(object):
    """
    The main object controlling the scene.

    This is the base Scene object.  It should be inherited, and then its 
    methods overridden.  
    """

    def __init__(self):
        """
        The init function
        """
        object.__init__(self)
        debugMsg("Called Scene.__init__()")

        self.renderer = Renderer()

        self.xSize = 640
        self.ySize = 480

    def add(self, obj):
        """
        Add a new item to the scene

        @param obj: The object to add to the scene
        @type obj: object
        """
        debugMsg("Called Scene.add()")

        if obj is None:
            raise ValueError, "You must specify an object to add"

        # print a warning message if get to here
        overrideWarning("Scene.add")

        return

    def delete(self, obj):
        """
        Delete an item from the scene

        @param obj: The object to remove
        @type obj: object
        """
        debugMsg("Called Scene.delete()")

        if obj is None:
            raise ValueError, "You must specify an object to delete"

        # print a warning message if get to here
        overrideWarning("Scene.delete")

        return

    def place(self, obj):
        """
        Place an object within a scene

        @param obj: The object to place within the scene
        @type obj: object
        """
        debugMsg("Called Scene.place()")

        if obj is None:
            raise ValueError, "You must specify an object to place"

        # print a warning message if get to here
        overrideWarning("Scene.place")

        return

    def render(self, pause=False, interactive=False):
        """
        Render (or re-render) the scene
        
        Render the scene, either to screen, or to a buffer waiting for a save

        @param pause: Flag to wait at end of script evaluation for user input
        @type pause: boolean
        
        @param interactive: Whether or not to have interactive use of the output
        @type interactive: boolean
        """
        debugMsg("Called Scene.render()")
        renderer = self.renderer

        # I don't yet know where to put this, but just to get stuff going...
        renderer.runString("# Scene.render()\n")

        # optionally print out the evaluation stack to make sure we're doing
        # the right thing
        debugMsg("Here is the evaluation stack")
        debugMsg(60*"#")
        debugMsg(renderer.getEvalStack())
        debugMsg(60*"#")

	# execute the eval stack
	evalStack = renderer.getEvalStack()
	exec evalStack in self.renderer.renderDict

        # flush the evaluation stack
        debugMsg("Flusing evaluation stack")
        renderer.resetEvalStack()

        # this is just to stop lint from complaining that pause and
        # interactive aren't used
        if pause is not True or pause is not False:
            raise ValueError, "\'pause\' must be either True or False"

        if interactive is not True or pause is not False:
            raise ValueError, "\'interactive\' must be either True or False"

        return

    def save(self, fname, format):
        """
        Save the scene to a file

        @param fname: The name of the file to save the scene to
        @type fname: string

        @param format: The format in which to save the scene
        @type format: Image object or string
        """
        debugMsg("Called Scene.save()")

        if fname is None or fname == "":
            raise ValueError, "You must specify an output filename"

        if format is None or format == "":
            raise ValueError, "You must specify an image format"

        # now check the type of arguments sent in
        if __debug__:
            assert isinstance(fname, str),  "Incorrect data type; expected str"
            assert isinstance(format, str), "Incorrect data type; expected str"

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("Scene.save")

        return

    write = save

    def setBackgroundColor(self, *color):
        """
        Sets the background color of the Scene

        @param color: The color to set the background to.  Can be RGB or CMYK
        @type color: tuple
        """
        debugMsg("Called Scene.setBackgroundColor()")

        # print a warning message if get to here
        overrideWarning("Scene.setBackgroundColor")

        # pity this code doesn't work....
        # need to check on the values given in the *color array.
        # if they're greater than 1, scale so that the largest is 1
        #maxColor = None
        #for i in range(len(color)):
            #if color[i] > 1:
                #maxColor = color[i]
                #print maxColor

        ## if a maximum colour is found, then scale the colours
        #if maxColor is not None:
            #for i in range(len(color)):
                #color[i] = color[i]/maxColor
        
        # if color is of length 3, then we have rgb
        # if length is 4 then cmyk
        # if length is 1 then greyscale
        # otherwise barf
        if len(color) == 3:
            # ok, using rgb
            # probably should use a Color object or something
            # this will do in the meantime
            pass
        else:
            raise ValueError, "Sorry, only RGB color is supported at present"

        return

    def getBackgroundColor(self):
        """
        Gets the current background color setting of the Scene
        """
        debugMsg("Called Scene.getBackgroundColor()")

        # print a warning message if get to here
        overrideWarning("Scene.getBackgroundColor")

        return

    def setSize(self, xSize, ySize):
        """
        Sets the size of the scene.

        This size is effectively the renderer window size.

        @param xSize: the size to set the x dimension
        @type xSize: int

        @param ySize: the size to set the y dimension
        @type ySize: int
        """
        debugMsg("Called Scene.setSize()")

        # make sure that the arguments are the right kind of thing
        if __debug__:
            assert isinstance(xSize, int), "Incorrect data type; expected int"
            assert isinstance(ySize, int), "Incorrect data type; expected int"

        self.xSize = xSize
        self.ySize = ySize

        return

    def getSize(self):
        """
        Gets the current size of the scene

        This size is effectively the renderer window size.  Returns a tuple
        of the x and y dimensions respectively, in pixel units(??).
        """
        debugMsg("Called Scene.getSize()")

        return (self.xSize, self.ySize)

    def rendererCommand(self, command):
        """
        Allows the user to run a low-level renderer-specific command directly

        @param command: The renderer command to run as a string
        @type command: string
        """
        debugMsg("Called Scene.rendererCommand()")
        # check that we get a string as input
        if __debug__:
            assert isinstance(command, str), "Incorrect data type; expected str"

        evalString = "%s" % command
        self.renderer.runString(evalString)
        return

# vim: expandtab shiftwidth=4:
