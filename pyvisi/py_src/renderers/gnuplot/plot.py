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

# $Id: plot.py,v 1.61 2005/11/30 04:28:30 paultcochrane Exp $

"""
Class and functions associated with a pyvisi Plot objects (gnuplot)
"""

# generic imports
from pyvisi.renderers.gnuplot.common import debugMsg
import copy

# module specific imports
from pyvisi.renderers.gnuplot.item import Item

__revision__ = '$Revision: 1.61 $'

class Plot(Item):
    """
    Abstract plot class
    """
    def __init__(self, scene):
        """
        Initialisation of abstract Plot class

        @param scene: the scene with which to associate the Plot
        @type scene: Scene object
        """
        debugMsg("Called Plot.__init__()")
        Item.__init__(self)

        self.renderer = scene.renderer

        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to Plot

        @param dataList: the data to set to the plot (should be an array or list
        or something)
        @type dataList: tuple
        """
        debugMsg("Called Plot.setData()")

        if dataList is None:
            raise ValueError, "You must specify a data list"
        
        return

    def setTitle(self, title):
        """
        Set the plot title

        @param title: the string holding the title to the plot
        @type title: string
        """
        debugMsg("Called Plot.setTitle()")

        self.title = title

        return

    def setXLabel(self, label):
        """
        Set the label of the x-axis

        @param label: the string holding the label of the x-axis
        @type label: string
        """
        debugMsg("Called Plot.setXLabel()")

        self.xlabel = label

        return

    def setYLabel(self, label):
        """
        Set the label of the y-axis

        @param label: the string holding the label of the y-axis
        @type label: string
        """
        debugMsg("Called Plot.setYLabel()")

        self.ylabel = label

        return

    def setZLabel(self, label):
        """
        Set the label of the z-axis

        @param label: the string holding the label of the z-axis
        @type label: string
        """
        debugMsg("Called Plot.setZLabel()")

        self.zlabel = label

        return

    def setLabel(self, axis, label):
        """
        Set the label of a given axis

        @param axis: string (Axis object maybe??) of the axis (e.g. x, y, z)

        @param label: string of the label to set for the axis
        @type label: string
        """
        debugMsg("Called Plot.setLabel()")

        # string-wise implementation (really budget implementation too)
        if axis == 'x' or axis == 'X':
            self.xlabel = label
        elif axis == 'y' or axis == 'Y':
            self.ylabel = label
        elif axis == 'z' or axis == 'Z':
            self.zlabel = label
        else:
            raise ValueError, "axis must be x or y or z"

        return

    def setXAxisRange(self, xMin, xMax):
        """
        Set the range of the x axis

        @param xMin: the minimum value of the x-axis
        @type xMin: float

        @param xMax: the maximum value of the x-axis
        @type xMax: float
        """
        debugMsg("Called Plot.setXAxisRange()")

        # check that the max is (strictly) bigger than the min
        if (xMax <= xMin):
            raise ValueError, "xMax should be strictly greater than xMin"

        evalString = "_gnuplot('set xrange [%f:%f]')" % (xMin, xMax)
        self.renderer.runString(evalString)

        return

    def setYAxisRange(self, yMin, yMax):
        """
        Set the range of the y axis

        @param yMin: the minimum value of the y-axis
        @type yMin: float

        @param yMax: the maximum value of the y-axis
        @type yMax: float
        """
        debugMsg("Called Plot.setYAxisRange()")

        # check that the max is (strictly) bigger than the min
        if (yMax <= yMin):
            raise ValueError, "yMax should be strictly greater than yMin"

        evalString = "_gnuplot('set yrange [%f:%f]')" % (yMin, yMax)
        self.renderer.runString(evalString)

        return

    def setZAxisRange(self, zMin, zMax):
        """
        Set the range of the z axis

        @param zMin: the minimum value of the z-axis
        @type zMin: float

        @param zMax: the maximum value of the z-axis
        @type zMax: float
        """
        debugMsg("Called Plot.setZAxisRange()")

        # check that the max is (strictly) bigger than the min
        if (zMax <= zMin):
            raise ValueError, "zMax should be strictly greater than zMin"

        evalString = "_gnuplot('set zrange [%f:%f])" % (zMin, zMax)
        self.renderer.runString(evalString)

        return

# vim: expandtab shiftwidth=4:

