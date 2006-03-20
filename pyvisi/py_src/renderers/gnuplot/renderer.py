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

# $Id: renderer.py,v 1.17 2005/11/07 04:42:08 paultcochrane Exp $

## @file render.py

"""
This is the file for the Renderer class
"""

# generic imports
from common import debugMsg
from common import _rendererName
from esys.pyvisi.renderer import Renderer as BaseRenderer

__revision__ = '$Revision'

class Renderer(BaseRenderer):
    """
    A generic object holding a renderer of a Scene().
    """

    def __init__(self):
        """
        Initialisation of Renderer() class
        """
        debugMsg("Called Renderer.__init__()")
        BaseRenderer.__init__(self)

        # initialise some attributes
        self.renderWindowWidth = 640
        self.renderWindowHeight = 480

        # what is the name of my renderer?
        self.name = _rendererName

        # initialise the evalstack
        self._evalStack = ""

        # the namespace to run the exec code
        self.renderDict = {}

        # keep the initial setup of the module for later reuse
        self._initStack = ""

        # initialise the renderer module
        self.addToInitStack("# Renderer.__init__()")
        self.addToInitStack("import Gnuplot")
        # Gnuplot package needs Numeric package
        self.addToInitStack("from Numeric import *")

        # need to add a check here to see if the Numeric module has been
        # imported, and if not, throw an error. 

        self.addToInitStack("_gnuplot = Gnuplot.Gnuplot()")

# vim: expandtab shiftwidth=4:
