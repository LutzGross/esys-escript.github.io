"""
Classes defining geometrical items in visualization

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


def Position(object):
    """
    A position in global coordinates
    """
    pass

def Origin(Position):
    """
    The position of the origin
    """
    pass

def Direction(object):
    """
    A dirction in global coordinates
    """
    pass

def XAxis(Direction):
    """
    The direction of the x-axis
    """
    pass

def YAxis(Direction):
    """
    The direction of the y-axis
    """
    pass

def ZAxis(Direction):
    """
    The direction of the z-axis
    """
    pass

def Plane(object):
    """
    A plane in global coordinates
    """
    pass

def XYPlane(Plane):
    """
    The XY plane orthogonal to the z-axis
    """
    pass

def YZPlane(Plane):
    """
    The YZ plane orthogonal to the x-axis
    """
    pass

def ZXPlane(Plane):
    """
    The ZX plane orthogonal to the y-axis
    """
    pass

def Sphere(object):
    """
    A sphere
    """
    pass

