"""
Class and functions associated with a pyvisi Scene

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


# generic imports
from common import debugMsg
from esys.pyvisi.scene import Scene as BaseScene

# module specific imports
from renderer import Renderer

import os, math

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

        self.xSize = 640
        self.ySize = 480

        # these variables are used to determin the location of the camera
        self.centre = None
        self.bounds = None

        self.objectList = []

    def add(self, obj):
        """
        Add a new item to the scene

        @param obj: The object to add to the scene
        @type obj: object
        """
        debugMsg("Called Scene.add()")

        if obj is None:
            raise ValueError, "You must specify an object to add"

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
            raise ValueError, "You must specify an object to place" 

        return

    def render(self, pause=False, interactive=False, save=False):
        """
        Render (or re-render) the scene
        
        Render the scene, either to screen, or to a buffer waiting for a save

        @param pause: Flag to wait at end of script evaluation for user input
        @type pause: boolean

        @param interactive: Whether or not to have interactive use of the output
        @type interactive: boolean
        """
        debugMsg("Called Scene.render()")
        self.__render(pause,interactive,False)

    def save(self, fname, format):
        """
        Save the scene to a file

        Possible formats are:
            - PNG
            - TGA (uncompressed targa)
            - CTGA (compressed targa)
            - PPM
            - SYS (e.g. BMP on Windows, PICT on MacOS)

        @param fname: the name of the file to save to
        @type fname: string

        @param format: the image format of the output file
        @type format: Image object or string
        """
        debugMsg("Called Scene.save()")
        self.__render(False,Flase,True,fname, format)
       
    def __render(self, pause=False, interactive=False, save=False, fname="out", format="ppm"):
        """
        Render (or re-render) the scene
        
        Render the scene, either to screen, or to a buffer waiting for a save

        @param pause: Flag to wait at end of script evaluation for user input
        @type pause: boolean

        @param interactive: Whether or not to have interactive use of the output
        @type interactive: boolean

        @param save: Whether or not to save the output when render
        @type save: boolean

        @param fname: the name of the file to save to
        @type fname: string

        @param format: the image format of the output file
        @type format: Image object or string
        """
        debugMsg("Called Scene.render()")
        renderer = self.renderer
        self.fname = fname
        # if the format is passed in as a string or object, react
        # appropriately
        if save:
           import types
           if type(format) is types.StringType:
               fmt = format.lower()
           else:
               fmt = format.format

           if fmt == "png":
               self.renderer.addToInitStack("Output_File_Type=N")
           elif fmt == "tga":
               self.renderer.addToInitStack("Output_File_Type=T")
           elif fmt == "ppm":
               self.renderer.addToInitStack("Output_File_Type=P")
           elif fmt == "ctga":
               self.renderer.addToInitStack("Output_File_Type=C")
           elif fmt == "sys":
               self.renderer.addToInitStack("Output_File_Type=S")
           else:
               raise ValueError, "Unknown graphics format.  I got %s" % fmt

        # this is the string from which to make the ini file
        renderer.addToInitStack("; PyVisi Povray renderer module ini file")

        if pause:
            renderer.addToInitStack("Pause_When_Done=on")
        else:
            renderer.addToInitStack("Pause_When_Done=off")

        if interactive:
            print "Interactive mode not available with POVRAY renderer module"

        # add the standard headers
        evalString = "#include \"shapes.inc\"\n"
        evalString += "#include \"colors.inc\"\n"

        # add a camera  (this should be its own object, and I should just
        # call some method of its own to get the info needed)

        # given the centre and bounds variables, determine the location and
        # look_at variables

        # do some checking first
        if self.centre is None or self.bounds is None:
            raise ValueError, "Model centre or bounds not defined."

        # calculate the y height
        yHeight = self.bounds[3] - self.bounds[2]
        xWidth = self.bounds[1] - self.bounds[0]
        if yHeight > xWidth:
            dist = yHeight
        else:
            dist = xWidth
        zDepth = self.bounds[5] - self.bounds[4]
        angle = 67.380*math.pi/180.0  # povray default
        opp = dist*1.5/2.0
        adj = opp/math.tan(angle/2.0)
        loc = adj + zDepth/2.0 + self.centre[2]

        # now write the camera to file
        evalString += "camera {\n"
        evalString += "  location <%f, %f, -%f>\n" % \
                (self.centre[0], self.centre[1], loc)
        evalString += "  direction <0, 0, 1>\n"
        evalString += "  up <0, 1, 0>\n"
        evalString += "  right <4/3, 0, 0>\n"
        evalString += "  look_at <%f, %f, -%f>\n" % \
                (self.centre[0], self.centre[1], self.centre[2])
        evalString += "}\n"

        # add the light source
        evalString += "light_source {\n"
        evalString += "  <0, 0, -500>\n"
        evalString += "  colour White\n"
        evalString += "}\n"

        renderer.runString(evalString)

        # write the resolution settings
        iniString = "Width=%d\n" % self.xSize
        iniString += "Height=%d\n" % self.ySize

        # might as well have antialiasing going
        iniString += "Antialias=on"

        renderer.addToInitStack(iniString)

        # get all objects in the scene to render themselves
        for obj in self.objectList:
            obj.render()

        # if saving to file don't render to the screen
        if save:
            renderer.addToInitStack("Output_To_File=on")
            renderer.addToInitStack("Display=off")
        else:
            renderer.addToInitStack("Output_To_File=off")
            renderer.addToInitStack("Display=on")

        # optionally print out the evaluation stack to make sure we're doing
        # the right thing
        debugMsg("Here is the evaluation stack")
        debugMsg(60*"#")
        debugMsg(renderer.getEvalStack())
        debugMsg(60*"#")

        if save:
            import re
            r = re.compile(r'\.\w+')
            self.fname = r.sub('', self.fname)
            povFname = "%s.pov" % self.fname
            iniFname = "%s.ini" % self.fname
        else:
            import time
            timeString = str(time.clock())
            povFname = "povPlot%s.pov" % timeString
            iniFname = "povPlot%s.ini" % timeString

        renderer.addToInitStack("Input_File_Name=%s" % povFname)

        ### generate the pov file
        pov = open(povFname, "w")
        pov.write(renderer.getEvalStack())
        pov.close()

        ### generate the ini file
        ini = open(iniFname, "w")
        ini.write(renderer.getInitStack())
        ini.close()

        # now compile the string object into code, and execute it
        if os.system("povray %s" % iniFname) != 0:
            print "evalStack execution failed"
            print "evalStack = \'%s\'" % renderer.getEvalStack()
            return None

        # flush the evaluation stack
        debugMsg("Flusing evaluation stack")
        renderer.resetEvalStack()

        # flush the init stack
        debugMsg("Flushing init stack")
        renderer.resetInitStack()

        # clean up a bit
        os.unlink(povFname)
        os.unlink(iniFname)

        return

    # set up an alias for the save method
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
        Gets the current background colour/color setting of the Scene
        """
        debugMsg("Called Scene.getBackgroundColor()")
        return

# vim: expandtab shiftwidth=4:
