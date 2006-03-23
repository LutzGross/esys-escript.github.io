"""
This is the file of the base class for items within a scene

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


from common import debugMsg, overrideWarning

class Item(object):
    """
    This is the base class for items within a scene
    """

    def __init__(self):
        """
        Initialisation
        """
        object.__init__(self)
        debugMsg("Called Item.__init__()")

    def render(self):
        """
        Render the object
        """
        debugMsg("Called Item.render()")

        # print a warning if get to here
        overrideWarning("Item.render")

        return

# vim: expandtab shiftwidth=4:
