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

# $Id: renderer.py,v 1.7 2005/11/07 04:53:34 paultcochrane Exp $

## @file renderer.py

"""
Brief introduction to what the file contains/does
"""

from pyvisi.renderers.plplot.common import debugMsg
from pyvisi.renderers.gnuplot.common import _rendererName
from pyvisi.renderer import Renderer as BaseRenderer

__revision__ = '$Revision: 1.7 $'

    
class Renderer(BaseRenderer):
    """
    Generic object holding a renderer of a Scene()
    """

    def __init__(self):
        """
        Initialisation
        """
        debugMsg("Called Renderer.__init__()")
        BaseRenderer.__init__(self)  # initialisation of base class
    
        # initialise some attributes
        self.renderWindowWidth = 640
        self.renderWindowHeight = 480

        self.name = _rendererName

        # the namespace to run the exec code
        self.renderDict = {}

        # initialise the evalstack
        self._evalStack = ""

        # keep the initial setup of the module for later reuse
        self._initStack = ""

        # initialise the renderer module
        self.addToInitStack("# Renderer.__init__()")
        self.addToInitStack("import plplot")

        # we now need the Numeric package so can handle arrays better
        self.addToInitStack("from Numeric import *")

# vim: expandtab shiftwidth=4:
