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

# $Id: scatter_plot_3d.py,v 1.1 2005/11/30 04:29:53 paultcochrane Exp $

"""
Class and functions associated with a pyvisi ScatterPlot3D objects (gnuplot)
"""

# generic imports
from common import debugMsg
import copy

# module specific imports
from plot import Plot

__revision__ = '$Revision: 1.1 $'

class ScatterPlot3D(Plot):
    """
    Scatter Plot in three dimensions.

    This is like a surface plot, but using crosses, or points as the data
    points.
    """

    def __init__(self, scene):
        """
        Initialisation of ScatterPlot3D class

        @param scene: the scene with which to associate the ScatterPlot3D
        @type scene: Scene object
        """
        debugMsg("Called ScatterPlot3D.__init__()")
        Plot.__init__(self, scene)

        # grab the renderer
        self.renderer = scene.renderer

        # set up some of the attributes
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None

        # now add the object to the scene
        scene.add(self)

    def setData(self, *dataList):
        """
        Sets the data to the given plot object.

        @param dataList: list of data objects to plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in ScatterPlot3D()")

        self.renderer.runString("# ScatterPlot3D.setData()")

        # for the moment, make sure that there are three arrays
        if len(dataList) != 3:
            raise ValueError, "Must have three arrays as input (at present)"

        # do some sanity checks on the data
        xData = dataList[0]
        yData = dataList[1]
        zData = dataList[2]

        if len(xData.shape) != 1:
            raise ValueError, "x data array is not of the correct shape: %s"\
                    % xData.shape

        if len(yData.shape) != 1:
            raise ValueError, "y data array is not of the correct shape: %s"\
                    % yData.shape

        if len(zData.shape) != 2:
            raise ValueError, "z data array is not of the correct shape: %s"\
                    % zData.shape

        # share around the data
        ## the x data
        self.renderer.renderDict['_x'] = copy.deepcopy(xData)

        ## the y data
        self.renderer.renderDict['_y'] = copy.deepcopy(yData)

        ## the z data
        self.renderer.renderDict['_z'] = copy.deepcopy(zData)

        self.renderer.runString(\
                "_data = Gnuplot.GridData(_z, _x, _y, binary=1)")

        return

    def render(self):
        """
        Does ScatterPlot3D object specific rendering stuff
        """
        debugMsg("Called ScatterPlot3D.render()")

        self.renderer.runString("# ScatterPlot3D.render()")
        evalString = "_gnuplot('set style data points')"
        self.renderer.runString(evalString)

        #self.renderer.runString("_gnuplot('set surface')")

        # if a title is set, put it here
        if self.title is not None:
            evalString = "_gnuplot.title(\'%s\')" % self.title
            self.renderer.runString(evalString)

        # if an xlabel is set, add it
        if self.xlabel is not None:
            evalString = "_gnuplot.xlabel(\'%s\')" % self.xlabel
            self.renderer.runString(evalString)

        # if a ylabel is set add it
        if self.ylabel is not None:
            evalString = "_gnuplot.ylabel(\'%s\')" % self.ylabel
            self.renderer.runString(evalString)

        # if a zlabel is set add it
        if self.zlabel is not None:
            evalString = "_gnuplot('set zlabel \\'%s\\'')" % self.zlabel
            self.renderer.runString(evalString)

        # set up the evalString to use for plotting
        evalString = "_gnuplot.splot(_data)"
        self.renderer.runString(evalString)

        return


# vim: expandtab shiftwidth=4:
