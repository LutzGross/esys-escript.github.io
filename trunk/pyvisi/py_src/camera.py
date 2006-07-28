"""
Class and functions associated with cameras and views

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
__author__="Paul Cochrane, L. Gross"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision:$"
__date__="$Date:$"

from common import PyvisiObject

class Camera(PyvisiObject):
    """
    A camera
    """
    pass

class FrontView(Camera):
    pass

class BackView(Camera):
    pass

class TopView(Camera):
    pass

class BottomView(Camera):
    pass

class LeftView(Camera):
    pass

class RightView(Camera):
    pass

class IsometricView(Camera):
    pass
