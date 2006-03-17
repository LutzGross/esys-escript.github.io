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

# $Id: __init__.py,v 1.10 2005/02/24 04:04:38 paultcochrane Exp $

## @file __init__.py

"""
Initialisation of the pyvisi base package
"""

from pyvisi.common import _pyvisiVersion, _pyvisiRevision
print "This is PyVisi version %s-%s" % (_pyvisiVersion, _pyvisiRevision)

__author__ = 'Paul Cochrane'
__version__ = _pyvisiVersion
__revision__ = _pyvisiRevision

from pyvisi.item import Item
from pyvisi.renderer import Renderer
from pyvisi.scene import Scene
from pyvisi.plot import Plot, \
        ArrowPlot, ContourPlot, LinePlot, SurfacePlot
from pyvisi.camera import Camera
from pyvisi.image import Image, \
        JpegImage, PdfImage, PngImage, PnmImage, TiffImage
from pyvisi.text import Text
from pyvisi.axes import Axes
from pyvisi.plane import Plane

# vim: expandtab shiftwidth=4:
