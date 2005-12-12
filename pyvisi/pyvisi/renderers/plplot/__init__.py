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

# $Id: __init__.py,v 1.4 2005/11/30 07:11:47 paultcochrane Exp $

## @file __init__.py

"""
Initialisation of the plplot renderer module
"""
    
from pyvisi.renderers.plplot.common \
        import _rendererName, _rendererVersion, _rendererRevision
print "This is the \"%s\" renderer module version %s-%s" % \
        (_rendererName, _rendererVersion, _rendererRevision)

__author__ = 'Paul Cochrane'
__version__ = _rendererVersion
__revision__ = _rendererRevision

from pyvisi.renderers.plplot.item import Item
from pyvisi.renderers.plplot.scene import Scene
from pyvisi.renderers.plplot.renderer import Renderer
from pyvisi.renderers.plplot.axes import Axes
from pyvisi.renderers.plplot.camera import Camera
from pyvisi.renderers.plplot.image import Image, \
        JpegImage, PbmImage, PdfImage, PngImage, PnmImage, PsImage, TiffImage
from pyvisi.renderers.plplot.plane import Plane

# plotting stuff
from pyvisi.renderers.plplot.plot import Plot
from pyvisi.renderers.plplot.arrow_plot import ArrowPlot
from pyvisi.renderers.plplot.contour_plot import ContourPlot
from pyvisi.renderers.plplot.line_plot import LinePlot
from pyvisi.renderers.plplot.surface_plot import SurfacePlot

from pyvisi.renderers.plplot.text import Text

# vim: expandtab shiftwidth=4:
