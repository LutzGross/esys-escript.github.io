# Copyright (C) 2004 Paul Cochrane
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

# $Id$

## @file plot.py

"""
@brief Class and functions associated with a pyvisi Plot objects
"""

from scene import Scene
from item import Item
from common import _debug

class Plot(Scene):
    """
    @brief Abstract plot class
    """
    def __init__(self):
        if _debug: print "\tCalled Plot.__init__()"
        pass

    def setData(self,data):
        if _debug: print "\tCalled setData() in Plot()"
        return True

class ArrowPlot(Plot):
    """
    @brief Arrow field plot
    """
    def __init__(self):
        if _debug: print "\tCalled ArrowPlot.__init__()"
        pass

    def setData(self,data):
        if _debug: print "\tCalled setData() in ArrowPlot()"
        return True

class ContourPlot(Plot):
    """
    @brief Contour plot
    """
    def __init__(self):
        if _debug: print "\tCalled ContourPlot.__init__()"
        pass

    def setData(self,data):
        if _debug: print "\tCalled setData() in ContourPlot()"
        return True

class LinePlot(Plot):
    """
    @brief Line plot
    """
    def __init__(self):
        if _debug: print "\tCalled LinePlot.__init__()"
        pass

    def setData(self,data):
        if _debug: print "\tCalled setData() in LinePlot()"
        return True

# vim: expandtab shiftwidth=4:

