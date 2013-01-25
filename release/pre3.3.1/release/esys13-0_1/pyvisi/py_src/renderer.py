"""
This is the file for the base Renderer class

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

from common import debugMsg

class Renderer(object):
    """
    A generic object holding a renderer of a Scene().
    """

    def __init__(self):
        """
        Initialisation of Renderer() class
        """
        object.__init__(self)
        debugMsg("Called Renderer.__init__()")

        # initialise some attributes
        self.renderWindowWidth = 640
        self.renderWindowHeight = 480

	# the namespace to run the exec code
	self.renderDict = {}

        # initialise the evalstack
        self._evalStack = ""

        # keep the initial setup of the module for later reuse
        self._initStack = ""

        # initialise the renderer module

    def setRenderWindowWidth(self, width):
        """
        Sets the render window width
        
        @param width: The width of the render window
        @type width: int
        """
        debugMsg("Called Renderer.setRenderWindowWidth()")
        # check that the argument makes sense
        if __debug__:
            assert isinstance(width, int), "Incorrect data type; expected int"

        self.renderWindowWidth = width
        return

    def setRenderWindowHeight(self, height):
        """
        Sets the render window height

        @param height: The height of the render window
        @type height: int
        """
        debugMsg("Called Renderer.setRenderWindowHeight()")
        # check that the argument makes sense
        if __debug__:
            assert isinstance(height, int), "Incorrect data type; expected int"

        self.renderWindowHeight = height
        return

    def getRenderWindowWidth(self):
        """
        Gets the render window width
        """
        debugMsg("Called Renderer.getRenderWindowWidth()")
        return self.renderWindowWidth

    def getRenderWindowHeight(self):
        """
        Gets the render window height
        """
        debugMsg("Called Renderer.getRenderWindowHeight()")
        return self.renderWindowHeight

    def setRenderWindowDimensions(self, width, height):
        """
        Sets the render window dimensions

        @param width: the width of the render window
        @type width: int

        @param height: the height of the render window
        @type height: int
        """
        debugMsg("Called Renderer.setRenderWindowDimensions()")
        # check that the argument makes sense
        if __debug__:
            assert isinstance(width, int), "Incorrect data type; expected int"
            assert isinstance(height, int), "Incorrect data type; expected int"

        self.renderWindowWidth = width
        self.renderWindowHeight = height

        return

    def getRenderWindowDimensions(self):
        """
        Gets the render window dimensions

        @return: tuple of window width and window height, respectively
        """
        debugMsg("Called Renderer.getRenderWindowDimensions()")
        return (self.renderWindowWidth, self.renderWindowHeight)

    def getEvalStack(self):
        """
        Gets the evaluation stack as it currently stands
        """
        debugMsg("Called Renderer.getEvalStack()")
        return self._evalStack

    def addToEvalStack(self, evalString):
        """
        Method to add commands to the evaluation stack
        
        @param evalString: The string of commands to be added to the evalStack
        @type evalString: string
        """
        debugMsg("Called Renderer.addToEvalStack()")
        # check that the argument is ok
        if __debug__:
            assert isinstance(evalString, str), \
                    "Incorrect data type; expected string"

        self._evalStack += evalString + '\n'
        return

    def runString(self, evalString):
        """
        Method to run the given string in the renderer python interpreter
        
        @param evalString: The string of commands to be run
        @type evalString: string
        """
        debugMsg("Called Renderer.runString()")
        # check that the argument is ok
        if __debug__:
            assert isinstance(evalString, str), \
                    "Incorrect data type; expected string"

        self._evalStack += evalString + '\n'
        return

    def getInitStack(self):
        """
        Gets the initialisation stack as it currently stands
        """
        debugMsg("called Renderer.getInitStack()")
        return self._initStack

    def addToInitStack(self, evalString):
        """
        Method to add commands to the reusable part of the evaluation stack

        @param evalString: The string of commands to be added to the evalStack
        @type evalString: string
        """
        debugMsg("Called Renderer.addToInitStack()")
        # check that the argument is ok
        if __debug__:
            assert isinstance(evalString, str)

        self._initStack += evalString + '\n'
        return

    def resetEvalStack(self):
        """
        Reset/flush the evaluation stack
        """
        debugMsg("Called Renderer.resetEvalStack()")
        self._evalStack = ""
        return

    def resetInitStack(self):
        """
        Reset/flush the initialisation stack
        """
        debugMsg("Called Renderer.resetInitStack()")
        self._initStack = ""
        return

# vim: expandtab shiftwidth=4:
