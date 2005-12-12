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

# $Id: surface_plot.py,v 1.1 2005/11/30 07:12:40 paultcochrane Exp $

"""
This file contains all of the classes for the various plotting objects.
"""

from pyvisi.renderers.plplot.common import debugMsg

from pyvisi.renderers.plplot.plot import Plot

__revision__ = '$Revision: 1.1 $'

 class SurfacePlot(Plot):
    """
    Surface plot
    """

    def __init__(self, scene):
        """
        Initialisation of SurfacePlot class

        @param scene: the scene with which to associate the SurfacePlot
        @type scene: Scene object
        """
        debugMsg("Called SurfacePlot.__init__()")
        Plot.__init__(self, scene)

        # grab the renderer
        self.renderer = scene.renderer

        # set up some of the attributes
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None

        # to show contours of the surface on the bottom of the axes, set
        # this variable to True
        self.contours = False

        # now add the object to the scene
        scene.add(self)

    def setData(self, *dataList):
        """
        Sets the data to the given plot object.

        @param dataList: list of data objects to plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in SurfacePlot()")

        self.renderer.runString("# SurfacePlot.setData()")

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

        # pass the data around
        ## the x data
        self.renderer.renderDict['_x'] = copy.deepcopy(xData)

        ## the y data
        self.renderer.renderDict['_y'] = copy.deepcopy(yData)

        ## the z data
        self.renderer.renderDict['_z'] = copy.deepcopy(zData)

        # determine the min and max of x, y and z in world coordinates
        evalString = "_xMin = min(_x)\n"
        evalString += "_xMax = max(_x)\n"
        evalString += "_yMin = min(_y)\n"
        evalString += "_yMax = max(_y)\n"
        evalString += "_zMin = min(_z.flat)\n"
        evalString += "_zMax = max(_z.flat)"
        self.renderer.runString(evalString)

        # min and max of x and y variables in normalised coordinates
        # values are those recommended in plplot documentation for surface
        # plots, and are hardcoded here for the meantime.
        evalString = "_xMin2D = -2.5\n"
        evalString += "_xMax2D = 2.5\n"
        evalString += "_yMin2D = -2.5\n"
        evalString += "_yMax2D = 4.0"
        self.renderer.runString(evalString)

        # sides of the box in normalised coordinates, again, these are those
        # recommended by plplot
        evalString = "_basex = 2.0\n"
        evalString += "_basey = 4.0\n"
        evalString += "_height = 3.0"    # possible name clash???
        self.renderer.runString(evalString)

        # the angle to view the box (this needs to be set up by a Camera()
        # object in the future)
        evalString = "_alt = 45.0\n"
        evalString += "_az = 30.0"
        self.renderer.runString(evalString)

        return

    def render(self):
        """
        Does SurfacePlot object specific rendering stuff
        """
        debugMsg("Called SurfacePlot.render()")

        self.renderer.runString("# SurfacePlot.render()")
        # initialise plplot 
        self.renderer.runString("plplot.plinit()")

        # set up the viewport for plotting
        evalString = "plplot.plenv(_xMin2D,_xMax2D,_yMin2D,_yMax2D, 0, -2)"
        self.renderer.runString(evalString)

        # set up the window
        evalString = "plplot.plw3d(_basex, _basey, _height, "
        evalString += "_xMin, _xMax, _yMin, _yMax, _zMin, _zMax, "
        evalString += "_alt, _az)"
        self.renderer.runString(evalString)

        # if a title is not set, set it to a null string
        # (this will help keep plplot happy)
        if self.title is not None:
            evalString = "plplot.plmtex(\"t\", 1.0, 0.5, 0.5, \"%s\")" % \
                    self.title
            self.renderer.runString(evalString)

        # if an xlabel is not set, set it to a null string
        if self.xlabel is None:
            self.xlabel = ""

        # if a ylabel is not set, set it to a null string
        if self.ylabel is None:
            self.ylabel = ""

        # if a zlabel is not set, set it to a null string
        if self.zlabel is None:
            self.zlabel = ""

        # put the labels (if any) on the graph.
        evalString = "plplot.plbox3(\"bnstu\", \"%s\", 0.0, 0, " % self.xlabel
        evalString += "\"bnstu\", \"%s\", 0.0, 0, " % self.ylabel
        evalString += "\"bcdmnstuv\", \"%s\", 0.0, 0)" % self.zlabel
        self.renderer.runString(evalString)

        # plot it!
        ### note: I can put other shading options into this call
        evalString = "plplot.plsurf3d(_x, _y, _z, 0, ())"
        self.renderer.runString(evalString)

        # finish stuff off
        self.renderer.runString("plplot.plend()")

        # if contours is true, set the relevant option
        if self.contours:
            pass
            # I need to implement this!!!  And put it further up in render()

        return
 
# vim: expandtab shiftwidth=4:
