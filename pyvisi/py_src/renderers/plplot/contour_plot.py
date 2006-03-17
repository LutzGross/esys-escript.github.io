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

# $Id: contour_plot.py,v 1.1 2005/11/30 07:12:40 paultcochrane Exp $

"""
This file contains all of the classes for the various plotting objects.
"""

from pyvisi.renderers.plplot.common import debugMsg

from pyvisi.renderers.plplot.plot import Plot

__revision__ = '$Revision: 1.1 $'

class ContourPlot(Plot):
    """
    Contour Plots
    """

    def __init__(self, scene):
        """
        Initialisation of the ContourPlot class

        @param arg: The scene with which to associate the plot
        @type arg: Scene object
        """
        debugMsg("Called ContourPlot.__init__()")
        Plot.__init__(self)  # initialisation of base class
        
        self.renderer = scene.renderer

        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None

        self.linestyle = None  # pyvisi-defined linestyle
        self._linestyle = None # renderer-specific linestyle
    
        # now add the object to the scene
        scene.add(self)

    def setData(self, *dataList, **options):
        """
        Sets the data to the given plot object.

        @param dataList: list of data objects to plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in ContourPlot()")

        self.renderer.runString("# ContourPlot.setData()")

        # this is a really dodgy way to get the data into the renderer
        # I really have to find a better, more elegant way to do this

        # for the moment, make sure that there are three arrays
        if len(dataList) != 3:
            raise ValueError, "Must have three arrays as input (at present)"

        # do some sanity checks on the data
        xData = dataList[0]
        yData = dataList[1]
        zData = dataList[2]

        if len(xData.shape) != 1:
            raise ValueError, "x data array is not of correct shape: %s" % \
                    xData.shape

        if len(yData.shape) != 1:
            raise ValueError, "y data array is not of correct shape: %s" % \
                    yData.shape

        if len(zData.shape) != 2:
            raise ValueError, "z data array is not of correct shape: %s" % \
                    zData.shape

        # share around the data
        ## the x data
        self.renderer.renderDict['_x'] = copy.deepcopy(xData)

        ## the y data
        self.renderer.renderDict['_y'] = copy.deepcopy(yData)

        ## the z data
        self.renderer.renderDict['_z'] = copy.deepcopy(zData)

        # determine the min and max of x, y and z
        evalString = "_xMin = min(x)\n"
        evalString += "_xMax = max(x)\n"
        evalString += "_yMin = min(y)\n"
        evalString += "_yMax = max(y)\n"
        evalString += "_zMin = \n"
        evalString += "_zMax = \n"

        ### please complete me!!

        self.renderer.runString(evalString)

        return

    def render(self):
        """
        Does ContourPlot object specific rendering stuff
        """
        debugMsg("Called ContourPlot.render()")

        self.renderer.runString("# ContourPlot.render()")
        self.renderer.runString("_gnuplot('set contour base')")
        self.renderer.runString("_gnuplot('set view 0, 0, 1, 1')")
        self.renderer.runString("_gnuplot('set nosurface')") # gnuplot 3.7

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

        self.renderer.runString("_gnuplot('set pm3d')")

        # set up the evalString to use for plotting
        evalString = "_gnuplot.splot(_data)"
        self.renderer.runString(evalString)

        return


# vim: expandtab shiftwidth=4:
