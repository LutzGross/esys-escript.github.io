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

# $Id: scene.py,v 1.30 2005/11/07 04:42:08 paultcochrane Exp $

## @file scene.py

"""
Class and functions associated with a pyvisi Scene (gnuplot)
"""

# generic imports
from common import debugMsg
from esys.pyvisi.scene import Scene as BaseScene

# module specific imports
from renderer import Renderer

__revision__ = '$Revision: 1.30 $'

class Scene(BaseScene):
    """
    The main object controlling the scene.
    
    Scene object methods and classes overriding the BaseScene class.
    """

    def __init__(self):
        """
        The init function
        """
        debugMsg("Called Scene.__init__()")
        BaseScene.__init__(self)

        self.renderer = Renderer()
        self.objectList = []

        self.xSize = 640
        self.ySize = 480

    def add(self, obj):
        """
        Add a new item to the scene

        @param obj: The object to add to the scene
        @type obj: object
        """
        debugMsg("Called Scene.add()")

        # make sure there is an object passed in
        if obj is None:
            raise ValueError, "You must specify an object to add"

        self.renderer.runString("# Scene.add()")
        self.objectList.append(obj)

        return

    def place(self, obj):
        """
        Place an object within a scene

        @param obj: The object to place within the scene
        @type obj: object
        """
        debugMsg("Called Scene.place()")

        if obj is None:
            raise ValueError, "You must specify an object to add"

        return

    def render(self, pause=False, interactive=False):
        """
        Render (or re-render) the scene
        
        Render the scene, either to screen, or to a buffer waiting for a save

        @param pause: Flag to wait at end of script evaluation for user input
        @type pause: boolean

        @param interactive: Whether or not to have interactive use of the 
        output (not available in all renderer modules)
        @type interactive: boolean
        """
        debugMsg("Called Scene.render()")
        renderer = self.renderer

        # gnuplot doesn't support interactive stuff (I think)
        if interactive:
            print "Gnuplot does not support scene interaction"
            print "Setting interactive to false"
            interactive = False

        renderer.runString("# Scene.render()")

        # get object added to the scene to render itself
        for obj in self.objectList:
            obj.render()
        
        # add some code to pause after rendering if asked to
        if pause:
            renderer.runString("raw_input(\"Press enter to continue\")")

        # prepend the init stack to the eval stack
        self.renderer._evalStack = self.renderer._initStack + \
                self.renderer._evalStack

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

        return

    def save(self, fname, format):
        """
        Save the scene to a file

        Possible formats are:
            - JPEG
            - Postscript
            - PNG
            - PBM

        @param fname: Name of output file
        @type fname: string

        @param format: Graphics format of output file
        @type format: Image object or string
        """
        debugMsg("Called Scene.save()")
        self.renderer.runString("# Scene.save()")

        # if the format is passed in as a string or object, react
        # appropriately
        import types
        if type(format) is types.StringType:
            fmt = format.lower()
        else:
            fmt = format.format

        # set the output format
        if fmt == "ps":
            self.renderer.runString(\
                    "_gnuplot('set terminal postscript color')")
        elif fmt == "png":
            evalString = "_gnuplot('set terminal png size %d,%d')" % \
                    (self.xSize, self.ySize)
            self.renderer.runString(evalString)
        elif fmt == "pbm":
            self.renderer.runString(\
                    "_gnuplot('set terminal pbm color')")
            # set the size of the output
            # have to scale properly to do this
            xScale = float(self.xSize)/640.0
            yScale = float(self.ySize)/480.0
            evalString = "_gnuplot('set size %f,%f')" % \
                    (xScale, yScale)
            self.renderer.runString(evalString)
        elif fmt == "jpeg" or fmt == "jpg":
            evalString = "_gnuplot('set terminal jpeg size %d,%d')" % \
                    (self.xSize, self.ySize)
            self.renderer.runString(evalString)
        else:
            raise ValueError, "Unknown graphics format.  I got: %s" % fmt

        # set the output filename
        evalString = "_gnuplot('set output \"%s\"')" % fname
        self.renderer.runString(evalString)

        # now render the whole shebang (again)
        self.render()

        return

    # alias the save() method to write()
    write = save

    def setBackgroundColor(self, *color):
        """
        Sets the background color of the Scene

        @param color: The color to set the background to.  Can be RGB or CMYK
        @type color: tuple
        """
        debugMsg("Called Scene.setBackgroundColor()")

        # pity this code doesn't work....
        # need to check on the values given in the *color array.
        # if they're greater than 1, scale so that the largest is 1
        #maxColor = None
        #for i in range(len(color)):
            #if color[i] > 1:
                #maxColor = color[i]
                #print maxColor
#
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
        return

# vim: expandtab shiftwidth=4:
