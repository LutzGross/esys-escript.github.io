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

# $Id: arrow_plot.py,v 1.1 2005/11/30 04:29:53 paultcochrane Exp $

"""
Class and functions associated with a pyvisi ArrowPlot objects (gnuplot)
"""

# generic imports
from common import debugMsg
import copy

# module specific imports
from plot import Plot

__revision__ = '$Revision: 1.1 $'

class ArrowPlot(Plot):
    """
    Arrow field plot
    """
    def __init__(self, scene):
        """
        Initialisation of ArrowPlot class

        @param scene: the scene with which to associate the ArrowPlot
        @type scene: Scene object
        """
        debugMsg("Called ArrowPlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer

        self.title = None
        self.xlabel = None
        self.ylabel = None

        # now add the object to the scene
        scene.add(self)

    def setData(self, *dataList):
        """
        Sets the data to the given plot object.

        @param dataList: list of data objects to plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in ArrowPlot()")
        
        self.renderer.runString("# ArrowPlot.setData()")

        # do some sanity checking on the data
        if len(dataList) != 4:
            raise ValueError, \
                    "Must have four vectors as input: x, y, dx, dy"

        for i in range(len(dataList)):
            if len(dataList[i].shape) != len(dataList[0].shape):
                raise ValueError, "All arrays must be of the same shape"

        for i in range(len(dataList)):
            if len(dataList[i].shape) != 1 and len(dataList[i].shape) != 2:
                errorString = \
                        "Can only handle 1D or 2D arrays: dim=%d" % \
                        len(dataList[i].shape)
                raise ValueError, errorString

        for i in range(len(dataList)):
            if len(dataList[0]) != len(dataList[i]):
                raise ValueError, "Input vectors must all be the same length"

        # if we have 2D arrays as input, we need to flatten them to plot the
        # data properly
        if len(dataList[0].shape) == 1:
            xData = dataList[0]
            yData = dataList[1]
            dxData = dataList[2]
            dyData = dataList[3]
        elif len(dataList[0].shape) == 2:
            xData = dataList[0].flat
            yData = dataList[1].flat
            dxData = dataList[2].flat
            dyData = dataList[3].flat
        else:
            raise ValueError, "Input vectors can only be 1D or 2D"

        # x data
        self.renderer.renderDict['_x'] = copy.deepcopy(xData)

        # y data
        self.renderer.renderDict['_y'] = copy.deepcopy(yData)

        # dx data
        self.renderer.renderDict['_dx'] = copy.deepcopy(dxData)

        # dy data
        self.renderer.renderDict['_dy'] = copy.deepcopy(dyData)

        # set up the data to plot
        evalString = \
                "_data = Gnuplot.Data(_x, _y, _dx, _dy, with=\'vectors\')"
        self.renderer.runString(evalString)

        return

    def render(self):
        """
        Does ArrowPlot object specific rendering stuff
        """
        debugMsg("Called ArrowPlot.render()")

        self.renderer.runString("# ArrowPlot.render()")

        # if a title is set, put it here
        if self.title is not None:
            evalString = "_gnuplot.title(\'%s\')" % self.title
            self.renderer.runString(evalString)

        # if an xlabel is set, add it
        if self.xlabel is not None:
            evalString = "_gnuplot.xlabel(\'%s\')" % self.xlabel
            self.renderer.runString(evalString)

        # if a ylabel is set, add it
        if self.ylabel is not None:
            evalString = "_gnuplot.ylabel(\'%s\')" % self.ylabel
            self.renderer.runString(evalString)

        # set up the evalString to use for plotting
        evalString = "_gnuplot.plot(_data)"
        self.renderer.runString(evalString)

        return


# vim: expandtab shiftwidth=4:
