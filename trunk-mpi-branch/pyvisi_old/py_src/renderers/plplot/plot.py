"""
This file contains all of the classes for the various plotting objects.

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


from common import debugMsg

from plot import Plot

class Plot(Item):
    """
    Abstract plot class
    """

    def __init__(self, scene):
        """
        Initialisation of abstract plot class

        @param scene: The scene with which to associate the plot
        @type scene: Scene object
        """
        debugMsg("Called Plot.__init__()")
        Item.__init__(self)  # initialisation of base class
    
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
        
# vim: expandtab shiftwidth=4:
