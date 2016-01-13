"""
Initialisation of the opendx renderer module

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

    
from common import _rendererName, _rendererVersion, _rendererRevision
print "This is the \"%s\" renderer module version %s-%s" % (_rendererName, _rendererVersion, _rendererRevision)

from item import Item
from scene import Scene
from renderer import Renderer
from axes import Axes
from camera import Camera
from image import Image, JpegImage, PbmImage, PdfImage, PngImage, PnmImage, PsImage, TiffImage
from plane import Plane
from plot import Plot, ArrowPlot, ContourPlot, LinePlot
from text import Text

# vim: expandtab shiftwidth=4:
