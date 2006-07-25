"""
Brief introduction to what the file contains/does

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

from item import Item

class Text(Item):
    """
    Brief introduction to what the class does
    """

    def __init__(self, arg):
        """
        Brief description of the init function

        @param arg: a description of the argument
        @type arg: the type of the argument
        """
        debugMsg("Called Text.__init__()")
        Item.__init__(self)  # initialisation of base class
    
# vim: expandtab shiftwidth=4:
