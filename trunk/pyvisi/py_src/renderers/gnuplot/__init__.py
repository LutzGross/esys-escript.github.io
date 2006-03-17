# Copyright (C) 2004-2005 Paul Cochrane
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# $Id: __init__.py,v 1.16 2005/11/30 04:28:30 paultcochrane Exp $

## @file __init__.py

"""
Initialisation of the gnuplot renderer module
"""

from pyvisi.renderers.gnuplot.common \
        import _rendererName, _rendererVersion, _rendererRevision
print "This is the \"%s\" renderer module version %s-%s" % \
        (_rendererName, _rendererVersion, _rendererRevision)

__author__ = 'Paul Cochrane'
__version__ = _rendererVersion
__revision__ = _rendererRevision

from pyvisi.renderers.gnuplot.item import Item
from pyvisi.renderers.gnuplot.renderer import Renderer
from pyvisi.renderers.gnuplot.scene import Scene

# plotting stuff
from pyvisi.renderers.gnuplot.plot import Plot
from pyvisi.renderers.gnuplot.arrow_plot import ArrowPlot
from pyvisi.renderers.gnuplot.contour_plot import ContourPlot
from pyvisi.renderers.gnuplot.line_plot import LinePlot
from pyvisi.renderers.gnuplot.mesh_plot import MeshPlot
from pyvisi.renderers.gnuplot.offset_plot import OffsetPlot
from pyvisi.renderers.gnuplot.scatter_plot import ScatterPlot
from pyvisi.renderers.gnuplot.scatter_plot_3d import ScatterPlot3D
from pyvisi.renderers.gnuplot.surface_plot import SurfacePlot


from pyvisi.renderers.gnuplot.camera import Camera
from pyvisi.renderers.gnuplot.image import Image, \
        JpegImage, PbmImage, PdfImage, PngImage, PnmImage, PsImage, TiffImage
from pyvisi.renderers.gnuplot.text import Text
from pyvisi.renderers.gnuplot.axes import Axes
from pyvisi.renderers.gnuplot.plane import Plane

# vim: expandtab shiftwidth=4:
