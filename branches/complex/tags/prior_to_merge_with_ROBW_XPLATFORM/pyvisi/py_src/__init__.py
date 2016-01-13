"""
Initialisation of the pyvisi base package

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

from item import Item
from renderer import Renderer
from scene import Scene
from plot import Plot, ArrowPlot, ContourPlot, LinePlot, SurfacePlot
from camera import Camera
from image import Image, JpegImage, PdfImage, PngImage, PnmImage, TiffImage
from text import Text
from axes import Axes
from plane import Plane

# vim: expandtab shiftwidth=4:
