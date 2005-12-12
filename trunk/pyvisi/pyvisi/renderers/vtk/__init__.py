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

# $Id: __init__.py,v 1.14 2005/11/30 03:05:37 paultcochrane Exp $

## @file __init__.py

"""
Initialisation of vtk renderer module
"""

from pyvisi.renderers.vtk.common \
        import _rendererName, _rendererVersion, _rendererRevision
print "This is the \"%s\" renderer module version %s-%s" % \
    (_rendererName, _rendererVersion, _rendererRevision)

__author__ = 'Paul Cochrane'
__version__ = _rendererVersion
__revision__ = _rendererRevision

from pyvisi.renderers.vtk.item import Item
from pyvisi.renderers.vtk.renderer import Renderer
from pyvisi.renderers.vtk.scene import Scene

# plotting stuff
from pyvisi.renderers.vtk.plot import Plot
from pyvisi.renderers.vtk.arrow_plot import ArrowPlot
from pyvisi.renderers.vtk.arrow_plot_3d import ArrowPlot3D
from pyvisi.renderers.vtk.ball_plot import BallPlot
from pyvisi.renderers.vtk.contour_plot import ContourPlot
from pyvisi.renderers.vtk.ellipsoid_plot import EllipsoidPlot
from pyvisi.renderers.vtk.isosurface_plot import IsosurfacePlot
from pyvisi.renderers.vtk.line_plot import LinePlot
from pyvisi.renderers.vtk.offset_plot import OffsetPlot
from pyvisi.renderers.vtk.surface_plot import SurfacePlot


from pyvisi.renderers.vtk.camera import Camera
from pyvisi.renderers.vtk.image import Image, \
        JpegImage, PdfImage, PngImage, PnmImage, PsImage, TiffImage
from pyvisi.renderers.vtk.text import Text
from pyvisi.renderers.vtk.axes import Axes
from pyvisi.renderers.vtk.plane import Plane
from pyvisi.renderers.vtk.box import Box, ClipBox

# vim: expandtab shiftwidth=4:
