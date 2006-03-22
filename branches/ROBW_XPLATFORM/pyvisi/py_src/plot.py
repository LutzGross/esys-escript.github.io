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

# $Id: plot.py,v 1.19 2005/11/02 04:52:07 paultcochrane Exp $

## @file plot.py

"""
Base class and functions associated with a pyvisi Plot objects
"""

# generic imports
from common import debugMsg, overrideWarning
from item import Item

__revision__ = '$Revision: 1.19 $'

class Plot(Item):
    """
    Abstract plot class

    This is the abstract base class of all Plot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of abstract Plot class

        @param scene: the scene with which to associate the Plot
        @type scene: Scene object
        """
        Item.__init__(self)
        debugMsg("Called Plot.__init__()")

        if scene is None:
            raise ValueError, "You must specify a scene object"

        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None

    def setData(self, *dataList):
        """
        Set data to Plot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in Plot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("Plot.setData")

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

        @param axis: string (Axis object maybe??) of the axis (e.g. x, y, z,)
        @type axis: string or Axis object

        @param label: string of the label to set for the axis
        @type label: string
        """
        debugMsg("Called Plot.setLabel()")

        # string-wise implementation
        if axis == 'x' or axis == 'X':
            self.xlabel = label
        elif axis == 'y' or axis == 'Y':
            self.ylabel = label
        elif axis == 'z' or axis == 'Z':
            self.zlabel = label
        else:
            raise ValueError, "axis must be x or y or z"

        return

class ArrowPlot(Plot):
    """
    Arrow field plot

    This is the abstract base class of all ArrowPlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of ArrowPlot class

        @param scene: the scene with which to associate the ArrowPlot
        @type scene: Scene object
        """
        Plot.__init__(self, scene)
        debugMsg("Called ArrowPlot.__init__()")

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to ArrowPlot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in ArrowPlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("ArrowPlot.setData")

        return

class BallPlot(Plot):
    """
    Ball plot

    This is the abstract base class of all BallPlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of BallPlot class

        @param scene: the scene with which to associate the BallPlot
        @type scene: Scene object
        """
        Plot.__init__(self, scene)
        debugMsg("Called BallPlot.__init__()")

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to BallPlot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in BallPlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("BallPlot.setData")

        return

class ContourPlot(Plot):
    """
    Contour plot

    This is the abstract base class of all ContourPlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of ContourPlot class

        @param scene: the scene with which to associate the ContourPlot
        @type scene: Scene object
        """
        Plot.__init__(self, scene)
        debugMsg("Called ContourPlot.__init__()")

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to ContourPlot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in ContourPlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("ContourPlot.setData")

        return

class EllipsoidPlot(Plot):
    """
    Ellipsoid plot

    This is the abstract base class of all EllipsoidPlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of EllipsoidPlot class

        @param scene: the scene with which to associate the EllipsoidPlot
        @type scene: Scene object
        """
        Plot.__init__(self, scene)
        debugMsg("Called EllipsoidPlot.__init__()")

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to EllipsoidPlot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in EllipsoidPlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("EllipsoidPlot.setData")

        return

class IsosurfacePlot(Plot):
    """
    Isosurface plot

    This is the abstract base class of all IsosurfacePlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of IsosurfacePlot class

        @param scene: the scene with which to associate the IsosurfacePlot
        @type scene: Scene object
        """
        Plot.__init__(self, scene)
        debugMsg("Called IsosurfacePlot.__init__()")

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to IsosurfacePlot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in IsosurfacePlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("IsosurfacePlot.setData")

        return

class LinePlot(Plot):
    """
    Line plot

    This is the abstract base class of all LinePlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of LinePlot class

        @param scene: the scene with which to associate the LinePlot
        @type scene: Scene object
        """
        Plot.__init__(self, scene)
        debugMsg("Called LinePlot.__init__()")

        self.renderer = scene.renderer

        if scene is None:
            raise ValueError, "You must specify a scene object"

        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None

        self.linestyle = None   # pyvisi-defined linestyle
        self._linestyle = None  # renderer-specific linestyle

    def setData(self, *dataList):
        """
        Sets the data to the given plot object.

        @param dataList: list of data objects to plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in LinePlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("LinePlot.setData")

        return

    def render(self):
        """
        Does LinePlot object specific (pre) rendering stuff
        """
        debugMsg("Called LinePlot.render()")

        # print a warning message if get to here
        overrideWarning("LinePlot.render")

        return

    def setLineStyle(self, linestyle):
        """
        Sets the linestyle of the LinePlot

        Linestyles may be either a word in the Gnuplot style, or a symbol 
        shortcut in the Matlab style.  Some of the options do not have a
        Matlab equivalent but do have a Gnuplot equivalent, or vice versa.

        What this method does, is take the linestyles possible as defined by
        PyVisi, and then does some conversion as best it can to get the
        relevant output from (in this case) gnuplot.  

        Possible linestyles are:
            1. lines ('-')
            2. points ('o')
            3. linespoints ('-o')
            4. dots ('.')
            5. dotted (':')
            6. dashes ('--')
            7. dotdashes ('-.')

        @param linestyle: the style to use for the lines
        @type linestyle: string
        """
        debugMsg("Called LinePlot.setLineStyle()")

        self.linestyle = linestyle

        # print a warning if get to here
        overrideWarning("LinePlot.setLineStyle")

        return

    def getLineStyle(self):
        """
        Gets the current linestyle of the LinePlot

        @return: the linestyle as a string
        """
        debugMsg("Called LinePlot.getLineStyle()")

        return self.linestyle

class ScatterPlot(Plot):
    """
    Scatter plot

    This is the abstract base class of all ScatterPlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """
    def __init__(self, scene):
        """
        Initialisation of ScatterPlot class

        @param scene: the scene with which to associate the ScatterPlot
        @type scene: Scene object
        """
        Plot.__init__(self, scene)
        debugMsg("Called ScatterPlot.__init__()")

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to ScatterPlot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called ScatterPlot.setData()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("ScatterPlot.setData")

        return

class ScatterPlot3D(Plot):
    """
    Three dimensional scatter plot

    This is the abstract base class of all ScatterPlot3D objects.  Renderer
    modules must inherit and override the methods defined here.
    """

    def __init__(self, scene):
        """
        Intialisation of ScatterPlot3D class

        @param scene: the scene with which to associate the ScatterPlot3D
        @type scene: Scene object
        """
        debugMsg("Called ScatterPlot3D.__init__()")
        Plot.__init__(self, scene)

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to a ScatterPlot3D

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called ScatterPlot3D.setData()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("ScatterPlot3D.setData")

        return

    def render(self):
        """
        Perform ScatterPlot3D specific rendering stuff
        """
        debugMsg("Called ScatterPlot3D.render()")

        # print a warning message if get to here
        overrideWarning("ScatterPlot3D.render")

        return

class SurfacePlot(Plot):
    """
    Surface plot

    This is the abstract base class of all SurfacePlot objects.  Renderer
    modules must inherit and override the methods defined here.
    """

    def __init__(self, scene):
        """
        Intialisation of SurfacePlot class

        @param scene: the scene with which to associate the SurfacePlot
        @type scene: Scene object
        """
        debugMsg("Called SurfacePlot.__init__()")
        Plot.__init__(self, scene)

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to a SurfacePlot

        @param dataList: the data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called SufracePlot.setData()")

        if dataList is None:
            raise ValueError, "You must specify a data list"

        # print a warning message if get to here
        overrideWarning("SurfacePlot.setData")

        return

    def render(self):
        """
        Perform SurfacePlot specific rendering stuff
        """
        debugMsg("Called SurfacePlot.render()")

        # print a warning message if get to here
        overrideWarning("SurfacePlot.render")

        return

# vim: expandtab shiftwidth=4:

