"""
The classes associated with Planes

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
from common import debugMsg, overrideWarning
from item import Item

class Plane(Item):
    """
    Generic class for Plane objects
    """

    def __init__(self, scene):
        """
        Initialisation of the Plane object

        @param scene: the scene object within which the plane is
        @type scene: Scene object
        """
        Item.__init__(self)
        debugMsg("Called Plane.__init__()")

        self.renderer = scene.renderer

    def mapImageToPlane(self, image):
        """
        Maps an Image object onto a Plane object

        @param image: the image object to be mapped
        @type image: Image object
        """
        debugMsg("Called Plane.mapImageToPlane()")

        if image is None:
            raise ValueError, "You must specify an image object"

        # print a warning message if get to here
        overrideWarning("Plane.mapImageToPlane")

        return

    def render(self):
        """
        Perform Plane object specific (pre)rendering tasks
        """
        debugMsg("Called Plane.mapImageToPlane()")

        # print a warning message if get to here
        overrideWarning("Plane.render")

        return
    
# vim: expandtab shiftwidth=4:
