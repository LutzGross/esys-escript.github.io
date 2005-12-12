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

# $Id: __init__.py,v 1.8 2005/06/24 00:43:28 paultcochrane Exp $

## @file __init__.py

"""
Initialisation of the povray renderer module
"""

from pyvisi.renderers.povray.common \
        import _rendererName, _rendererVersion, _rendererRevision
print "This is the \"%s\" renderer module version %s-%s" % \
        (_rendererName, _rendererVersion, _rendererRevision)

__author__ = 'Paul Cochrane'
__version__ = _rendererVersion
__revision__ = _rendererRevision

from pyvisi.renderers.povray.item import Item
from pyvisi.renderers.povray.renderer import Renderer
from pyvisi.renderers.povray.scene import Scene
from pyvisi.renderers.povray.plot import Plot, \
        ArrowPlot, ArrowPlot3D, BallPlot, ContourPlot, LinePlot
from pyvisi.renderers.povray.camera import Camera
from pyvisi.renderers.povray.image import Image, \
        JpegImage
from pyvisi.renderers.povray.text import Text
from pyvisi.renderers.povray.axes import Axes
from pyvisi.renderers.povray.plane import Plane

# vim: expandtab shiftwidth=4:
