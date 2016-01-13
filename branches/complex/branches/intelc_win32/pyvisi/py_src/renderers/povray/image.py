"""
Class and functions associated with a pyvisi Image objects

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
from esys.pyvisi.common import fileCheck

# module specific imports
from item import Item

class Image(Item):
    """
    Image class.  Generic class to handle image data.
    """
    def __init__(self, format, scene):
        """
        Initialises the Image class object
        
        @param format: The image format
        @type format: string

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        debugMsg("Called Image.__init__()")
        Item.__init__()

        #if format == "jpeg":
            #if _debug: print "\t%s: Using jpeg image format" % rendererName
            #return JpegImage(scene)
        #else:
            #print "Unknown image format %s" % format
            #return None

        if format is None:
            raise ValueError, "You must specify an image format"

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def load(self, fname):
        """
        Loads image data from file.

        @param fname: The filename from which to load image data
        @type fname: string
        """
        debugMsg("Called Image.load()")

        fileCheck(fname)

        return

class JpegImage(Image):
    """
    Subclass of Image class to explicitly handle jpeg images
    """
    def __init__(self, scene):
        """
        Initialises the JpegImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        debugMsg("Called JpegImage.__init__()")
        Image.__init__()
        self.renderer = scene.renderer

    def load(self, fname):
        """
        Loads jpeg image data from file.

        @param fname: The filename from which to load jpeg image data
        @type fname: string
        """
        debugMsg("Called JpegImage.load()")

        fileCheck(fname)
        
        return

    def render(self):
        """
        Does JpegImage object specific (pre)rendering stuff
        """
        debugMsg("Called JpegImage.render()")

        return


# vim: expandtab shiftwidth=4:
