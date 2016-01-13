"""
Class and functions associated with a pyvisi Text object

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

class Text(Item):
    """
    Text
    """
    def __init__(self, scene):
        """
        Initialisation of the Text object

        @param scene: the scene with which to associate the Text object
        @type scene: Scene object
        """
        debugMsg("Called Text.__init__()")
        Item.__init__(self)
        self.font = "Times"

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setFont(self, font):
        """
        Set the current font

        @param font: the font to set
        @type font: string
        """
        self.font = font
        return

    def getFont(self):
        """
        Get the current font
        """
        return self.font

# vim: expandtab shiftwidth=4:
