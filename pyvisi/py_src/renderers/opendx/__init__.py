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

# $Id: __init__.py,v 1.2 2005/02/08 07:31:34 paultcochrane Exp $

## @file __init__.py

"""
Initialisation of the opendx renderer module
"""
    
from pyvisi.renderers.opendx.common \
        import _rendererName, _rendererVersion, _rendererRevision
print "This is the \"%s\" renderer module version %s-%s" % \
        (_rendererName, _rendererVersion, _rendererRevision)

__author__ = 'Paul Cochrane'
__version__ = _rendererVersion
__revision__ = _rendererRevision

from pyvisi.renderers.opendx.item import Item
from pyvisi.renderers.opendx.scene import Scene
from pyvisi.renderers.opendx.renderer import Renderer
from pyvisi.renderers.opendx.axes import Axes
from pyvisi.renderers.opendx.camera import Camera
from pyvisi.renderers.opendx.image import Image, \
        JpegImage, PbmImage, PdfImage, PngImage, PnmImage, PsImage, TiffImage
from pyvisi.renderers.opendx.plane import Plane
from pyvisi.renderers.opendx.plot import Plot, \
        ArrowPlot, ContourPlot, LinePlot
from pyvisi.renderers.opendx.text import Text

# vim: expandtab shiftwidth=4:
