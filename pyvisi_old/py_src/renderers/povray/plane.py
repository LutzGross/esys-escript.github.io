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
from common import debugMsg

# module specific imports
from item import Item

class Plane(Item):
    """
    Generic class for Plane objects
    """

    def __init__(self, scene):
        """
        Initialisation of the Plane object
        """
        debugMsg("Called Plane.__init__()")
        Item.__init__()

        self.renderer = scene.renderer

    def mapImageToPlane(self, image):
        """
        Maps an Image object onto a Plane object
        """
        debugMsg("Called Plane.mapImageToPlane()")

        if image is None:
            raise ValueError, "You must specify an image object"

        return

    def render(self):
        """
        Perform Plane object specific (pre)rendering tasks
        """
        debugMsg("Called Plane.mapImageToPlane()")

        return
    
# vim: expandtab shiftwidth=4:
