"""
Initialisation of the gnuplot renderer module

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
print "This is the \"%s\" renderer module version %s-%s"%(_rendererName, _rendererVersion, _rendererRevision)

__revision__ = _rendererRevision

from item import Item
from renderer import Renderer
from scene import Scene

# plotting stuff
from plot import Plot
from arrow_plot import ArrowPlot
from contour_plot import ContourPlot
from line_plot import LinePlot
from mesh_plot import MeshPlot
from offset_plot import OffsetPlot
from scatter_plot import ScatterPlot
from scatter_plot_3d import ScatterPlot3D
from surface_plot import SurfacePlot


from camera import Camera
from image import Image, JpegImage, PbmImage, PdfImage, PngImage, PnmImage, PsImage, TiffImage
from text import Text
from axes import Axes
from plane import Plane

# vim: expandtab shiftwidth=4:
