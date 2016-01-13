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

# $Id: plot.py,v 1.9 2005/11/02 04:52:08 paultcochrane Exp $

## @file plot.py

"""
Class and functions associated with a pyvisi Plot objects
"""

# generic imports
from common import debugMsg

# module specific imports
from item import Item

import Numeric
import math
import os

__revision__ = '$Revision: 1.9 $'

class Plot(Item):
    """
    Abstract plot class
    """
    def __init__(self, scene):
        """
        Initialisation of Plot class

        @param scene: the scene with which to associate the plot
        @type scene: Scene object
        """
        debugMsg("Called Plot.__init__()")
        Item.__init__(self)

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set the data to the plot

        @param dataList: list of data objects to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in Plot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"
        
        return

class ArrowPlot(Plot):
    """
    Arrow field plot
    """
    def __init__(self, scene):
        """
        Initialisation of the ArrowPlot class

        @param scene: the scene with which to associate the arrow plot
        @type scene: Scene object
        """
        debugMsg("Called ArrowPlot.__init__()")
        Plot.__init__()

        self.renderer = scene.renderer

    def setData(self, *dataList):
        """
        Set data to the plot

        @param dataList: list/tuple of data objects to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in ArrowPlot()")
        self.renderer.runString("// ArrowPlot.setData()")

        if dataList is None:
            raise ValueError, "You must specify a data list"
        
        return

    def render(self):
        """
        Does ArrowPlot specific rendering tasks
        """
        debugMsg("Called render() in ArrowPlot")
        self.renderer.runString("// ArrowPlot.render()")

        return

class ArrowPlot3D(Plot):
    """
    Arrow field plot in three dimensions
    """
    def __init__(self, scene):
        """
        Initialisation of the ArrowPlot3D class

        @param scene: the scene with which to associate the arrow plot
        @type scene: Scene object
        """
        debugMsg("Called ArrowPlot3D.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer
        self.renderer.runString("// ArrowPlot3D.__init__()")

        self.scene = scene

        # options
        self.fname = None
        self.format = None

        self.title = None

        # the data to pass between setData() and render()
        self.x = None
        self.y = None
        self.z = None
        self.dx = None
        self.dy = None
        self.dz = None

        self.norm = None
        self.maxNorm = None

        # share the colours around
        self.red = None
        self.green = None
        self.blue = None
 
        # add the plot to the scene
        scene.add(self)

    def setData(self, *dataList, **options):
        """
        Set data to the plot

        @param dataList: list/tuple of data objects to set to the plot
        @type dataList: tuple

        @param options: dictionary of extra options to method
        @type options: dict
        """
        debugMsg("Called setData() in ArrowPlot3D()")
        self.renderer.runString("// ArrowPlot3D.setData()")

        # process the options if any
        ## fname
        if options.has_key('fname'):
            fname = options['fname']
        else:
            fname = None
        ## format
        if options.has_key('format'):
            format = options['format']
        else:
            format = None

        # I think we want to pass this info around
        self.fname = fname
        self.format = format

        # do some sanity checking on the inputs
        if fname is None and format is not None:
            raise ValueError, "Format specified, but no input filename"
        elif fname is not None and format is None:
            raise ValueError, "Filename specified, but no format"
        elif (fname is not None or format is not None) and len(dataList) != 0:
            raise ValueError, \
                "Cannot specify a data list and an input file simultaneously"

        # of, if we have a data list and no args, use the data list
        if len(dataList) != 0 and fname is None and format is None:
            # do some sanity checking on the data
            if len(dataList) != 6:
                raise ValueError, \
                    "Must have six vectors as input: x,y,z,dx,dy,dz, found: %d"\
                            % len(dataList)
    
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
                    raise ValueError, \
                            "Input vectors must all be the same length"
    
            # if we have 2D arrays as input, we need to flatten them to plot the
            # data properly
            if len(dataList[0].shape) == 1:
                x = dataList[0]
                y = dataList[1]
                z = dataList[2]
                dx = dataList[3]
                dy = dataList[4]
                dz = dataList[5]
            elif len(dataList[0].shape) == 2:
                x = dataList[0].flat
                y = dataList[1].flat
                z = dataList[2].flat
                dx = dataList[3].flat
                dy = dataList[4].flat
                dz = dataList[5].flat
            else:
                raise ValueError, "Input vectors can only be 1D or 2D"

            # keep the number of points for future use
            numPoints = len(x)
    
            # calculate the norm of the data
            norm = Numeric.zeros(numPoints, typecode=Numeric.Float)
            for j in range(numPoints):
                norm[j] = math.sqrt(dx[j]*dx[j] + dy[j]*dy[j] + dz[j]*dz[j])

            # calculate the maximum norm of the data
            maxNorm = max(norm)

            # work out the centre of the model
            xCentre = (max(x) - min(x))/2.0
            yCentre = (max(y) - min(y))/2.0
            zCentre = (max(z) - min(z))/2.0
            centre = Numeric.array([xCentre, yCentre, zCentre])

            # work out the bounds of the model
            # (this should hopefully help set up the default camera)
            bounds = Numeric.array([min(x), max(x), min(y), max(y), \
                    min(z), max(z)])
 
        elif len(dataList) == 0 and fname is not None and format is not None:
            # well, lets process the vtk file then

            # had best make sure it exists
            if not os.path.exists(fname):
                raise SystemError, "File: %s doesn't exist" % fname

            # import vtk
            import vtk

            if format == 'vtk':
                # read old-style vtk files
                reader = vtk.vtkUnstructuredGridReader()
            elif format == 'vtk-xml':
                # read vtk xml files
                reader = vtk.vtkXMLUnstructuredGridReader()
            else:
                # barf
                raise ValueError, "Unknown format.  I got %s" % format

            reader.SetFileName(fname)
            reader.Update()

            # grab the grid
            grid = reader.GetOutput()

            # get the centre and bounds
            centre = grid.GetCenter()
            bounds = grid.GetBounds()

            # grab the point data
            points = grid.GetPoints()
            numPoints = points.GetNumberOfPoints()
            x = Numeric.zeros(numPoints, typecode=Numeric.Float)
            y = Numeric.zeros(numPoints, typecode=Numeric.Float)
            z = Numeric.zeros(numPoints, typecode=Numeric.Float)
            for i in range(numPoints):
                x[i], y[i], z[i] = points.GetPoint(i)

            # get the data at the points
            data = grid.GetPointData().GetVectors()
            dx = Numeric.zeros(numPoints, typecode=Numeric.Float)
            dy = Numeric.zeros(numPoints, typecode=Numeric.Float)
            dz = Numeric.zeros(numPoints, typecode=Numeric.Float)
            norm = Numeric.zeros(numPoints, typecode=Numeric.Float)
            for i in range(numPoints):
                dx[i], dy[i], dz[i] = data.GetTuple3(i)
                norm[i] = math.sqrt(dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i])

            # get the maximum norm
            maxNorm = grid.GetPointData().GetVectors().GetMaxNorm()

            # make a lookup table and invert it
            lut = vtk.vtkLookupTable()
            refLut = vtk.vtkLookupTable()
            lut.Build()
            refLut.Build()
            for j in range(256):
                lut.SetTableValue(j, refLut.GetTableValue(255-j))

            # get the colours
            red = Numeric.zeros(numPoints, typecode=Numeric.Float)
            green = Numeric.zeros(numPoints, typecode=Numeric.Float)
            blue = Numeric.zeros(numPoints, typecode=Numeric.Float)
            for i in range(numPoints):
                red[i], green[i], blue[i] = lut.GetColor(norm[i]/maxNorm)

            # share the colours around
            ### note: will need to work out how to do colourmaps myself 
            ### one day so that povray can have colours and just have data
            ### piped into it
            self.red = red
            self.green = green
            self.blue = blue
 
        else:
            # barf
            print len(dataList)
            print fname
            print format
            raise ValueError, "Shouldn't have got to here."

        self.scene.centre = centre
        self.scene.bounds = bounds

        # share the data around
        self.x = x
        self.y = y
        self.z = z
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.norm = norm
        self.maxNorm = maxNorm

        return

    def render(self):
        """
        Does ArrowPlot3D specific rendering tasks
        """
        debugMsg("Called render() in ArrowPlot3D")
        self.renderer.runString("// ArrowPlot3D.render()")

        # grab the data
        x = self.x
        y = self.y
        z = self.z
        dx = self.dx
        dy = self.dy
        dz = self.dz
        norm = self.norm
        maxNorm = self.maxNorm

        # grab the colours
        red = self.red
        green = self.green
        blue = self.blue

        # declare the arrow object
        evalString = "# declare Arrow = union {\n"
        evalString += "  cone {\n"
        evalString += "    <0, 0, 0>, 0.3\n"
        evalString += "    <1, 0, 0>, 0.0\n"
        evalString += "  }\n"
        evalString += "  cylinder {\n"
        evalString += "    <-1, 0, 0>\n"
        evalString += "    <0, 0, 0>,\n"
        evalString += "    0.15\n"
        evalString += "  }\n"
        evalString += "}\n"

        self.renderer.runString(evalString)

        evalString = ""
        for i in range(len(x)):
            evalString += "object {\n"
            scale = 0.05*norm[i]/maxNorm
            if scale < 1e-8:
                scale = 1e-8
            if norm[i] < 1e-8:
                norm[i] = 1e-8
            evalString += "  Arrow scale %g " % scale
            # note that these are the WRONG angle transformations!!
            xAngle = math.acos(dx[i]/norm[i])*180.0/math.pi
            yAngle = math.acos(dy[i]/norm[i])*180.0/math.pi
            zAngle = math.acos(dz[i]/norm[i])*180.0/math.pi
            evalString += "rotate <%f, %f, %f> " % (xAngle, yAngle, zAngle)
            evalString += "translate <%f, %f, -%f> " % (x[i], y[i], z[i])
            if self.fname is not None and self.format is not None:
                evalString += "pigment { colour <%f, %f, %f> }\n" % \
                        (red[i], green[i], blue[i])
            else:
                evalString += "pigment { colour Red }\n"
            evalString += "}\n"

        self.renderer.runString(evalString)

        # set the title if set
        if self.title is not None:
            # not implemented yet
            pass

        return 

class BallPlot(Plot):
    """
    Ball plot
    """
    def __init__(self, scene):
        debugMsg("Called BallPlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer
        self.renderer.runString("// BallPlot.__init__()")

        self.scene = scene

        # the data to pass between setData() and render()
        self.x = None
        self.y = None
        self.z = None
        self.radii = None

        self.title = None

        # the colour info to pass between setData() and render()
        self.red = None
        self.green = None
        self.blue = None

        # add the plot to the scene
        scene.add(self)

    def setData(self, points=None, radii=None, tags=None, colors=None,
            fname=None, format=None):
        """
        Set data to the plot
        @param points: the array to use for the points of the sphere
        locations in space
        @type points: float array

        @param fname: the name of the input vtk file
        @type fname: string

        @param format: the format of the input vtk file ('vtk' or 'vtk-xml')
        @type format: string 

        @param radii: the name of the scalar array in the vtk unstructured
        grid to use as the radii of the balls
        @type radii: float array

        @param colors: the name of the scalar array in the vtk unstructured
        grid to use as the colour tags of the balls
        @type colors: string

        @param tags: the name of the scalar array in the vtk unstructured
        grid to use as the colour of the tags of the balls
        @type tags: integer array
        """
        debugMsg("Called setData() in BallPlot()")
        self.renderer.runString("// BallPlot.setData()")

        # check that we have enough info to start with
        if points is None and fname is None:
            errorString = "You must supply either appropriate arrays of\n\
                    data or the name of a vtk file for input"
            raise ValueError, errorString

        # need to also check if they're both specified
        if points is not None and fname is not None:
            raise ValueError, \
                    "Sorry, you can't specify both a data list and a filename"

        # now check the bits required if the fname option is set
        if fname is not None:
            if format is None:
                raise ValueError, "You must specify a vtk file format"
            elif radii is None:
                raise ValueError, \
                "You must specify the name of the scalars to use as the radius"
            elif colors is None and tags is None:
                raise ValueError, \
                "You must specify the name of the scalars to use for the colors"

        # now check that the format is logical
        if format is not None and format != "vtk" and format != "vtk-xml":
            errorString = \
                "Unknown format: must be either 'vtk' or 'vtk-xml'\n"
            errorString += "I got: %s" % format
            raise ValueError, errorString

        # for now just hope that if the stuff is specified, that it agrees
        # with what's in the vtk unstructured grid

        if format is not None and (format == "vtk-xml" or format == "vtk"):
            # we're using vtk to read in the data, so import the module
            import vtk

            if format == "vtk-xml":
                debugMsg("Using vtk-xml file as input")
                # create the reader of the file
                reader = vtk.vtkXMLUnstructuredGridReader()
                reader.SetFileName(fname)
                reader.Update()
            elif format == "vtk":
                debugMsg("Using old-style vtk file as input")
                # create the reader of the file
                reader = vtk.vtkUnstructuredGridReader()
                reader.SetFileName(fname)
                reader.Update()

            # read the output to an unstructured grid
            grid = reader.GetOutput()

            # get the centre and bounds
            centre = grid.GetCenter()
            bounds = grid.GetBounds()

            # grab the points where the data sit
            vtkPoints = grid.GetPoints()

            # note that these next few steps are only necessary in vtk 4.2,
            # 4.4 grab the data to use for the radii of the balls
            vtkRadii = grid.GetPointData().GetScalars(radii)
    
            # grab the data to use for colouring the balls
            vtkTags = grid.GetPointData().GetScalars(tags)
    
            # now work out the number of tags, and their values
            numPoints = vtkTags.GetNumberOfTuples()
            valueDict = {}
            for i in range(numPoints):
                tagValue = vtkTags.GetValue(i)
                valueDict[tagValue] = 1
    
            numTags = len(valueDict.keys())
    
            tagValues = valueDict.keys()
            tagValues.sort()
    
            # now count the number of tags, and make an evenly spaced
            # array of points between zero and one, then use these as the
            # scalars to colour by
            vtkScaledTags = vtk.vtkFloatArray()
            vtkScaledTags.SetNumberOfTuples(numPoints)
            vtkScaledTags.SetNumberOfComponents(1)
            vtkScaledTags.SetName("scaledTags")
            if numTags == 1:
                for i in range(numPoints):
                    vtkScaledTags.InsertTuple1(i, 0.0)
            else:
                for i in range(numPoints):
                    tagValue = vtkTags.GetValue(i)
                    for j in range(numTags):
                        if tagValues[j] == tagValue:
                            vtkScaledTags.InsertTuple1(i, \
                                    float(j)/float(numTags-1))

            # use vtk to generate the colour map, will have to do this
            # myself for the non-vtk loading version
            lut = vtk.vtkLookupTable()
            lut.Build()

            red = Numeric.zeros(numPoints, typecode=Numeric.Float)
            green = Numeric.zeros(numPoints, typecode=Numeric.Float)
            blue = Numeric.zeros(numPoints, typecode=Numeric.Float)
            for i in range(numPoints):
                red[i], green[i], blue[i] = \
                        lut.GetColor(vtkScaledTags.GetValue(i))

            # now convert the information we want (radii, colours,
            # positions) into array objects so that I can play with them as
            # per normal in python

            x = Numeric.zeros(numPoints, typecode=Numeric.Float)
            y = Numeric.zeros(numPoints, typecode=Numeric.Float)
            z = Numeric.zeros(numPoints, typecode=Numeric.Float)
            radii = Numeric.zeros(numPoints, typecode=Numeric.Float)
            scaledTags = Numeric.zeros(numPoints, typecode=Numeric.Float)
            tags = Numeric.zeros(numPoints, typecode=Numeric.Float)
            for i in range(numPoints):
                ### the points
                x[i], y[i], z[i] = vtkPoints.GetPoint(i)

                ### the radii
                radii[i] = vtkRadii.GetValue(i)

                ### the tags
                scaledTags[i] = vtkScaledTags.GetValue(i)
                tags[i] = vtkTags.GetValue(i)
    
        elif format is None and points is not None:
            debugMsg("Using user-defined point data in BallPlot.setData()")
            ### if we get to here, then we have to construct the points,
            ### the radii and the colours all from scratch, add them to the
            ### grid and get everthing in the same form that we have for the
            ### case where we load a vtk data file

            numPoints = len(points)
            numRadii = len(radii)
            # do some sanity checking on the data
            if numRadii != numPoints:
                raise ValueError, \
                    "The number of points does not equal the number of radii"

            ### get the x, y, z data from the points
            x = Numeric.zeros(numPoints, typecode=Numeric.Float)
            y = Numeric.zeros(numPoints, typecode=Numeric.Float)
            z = Numeric.zeros(numPoints, typecode=Numeric.Float)

            for i in range(numPoints):
                x[i], y[i], z[i] = points[i]

            # find the bounds and the centre of the model
            centre = Numeric.array([(max(x)-min(x))/2.0, \
                    (max(y)-min(y))/2.0, (max(z)-min(z))/2.0])
            bounds = Numeric.array([min(x), max(x), min(y), max(y), \
                    min(z), max(z)])

            # make the colours
            if colors is None:
                # ok then, since we have no info, make them related to the
                # radii
                pass

            # what if we have a tags argument as well, and use that for the
            # colours if the colours array doesn't exists (which would be
            # more complicated for the user to set up, but can have the
            # functionality if someone wants to use it)

            if tags is None:
                debugMsg("Autogenerating tags in BallPlot.setData()")
                # relate the tags to the radii
                # need to find the number of different radii
                radiiDict = {}
                for i in range(numPoints):
                    radiiDict[str(radii[i])] = 1
                numRadii = len(radiiDict.keys())
                radiiKeys = radiiDict.keys()
                # now just make a list of evenly spaced tags up to numRadii
                tagValues = range(numRadii)
                numTags = numRadii

                tags = Numeric.zeros(numPoints, typecode=Numeric.Int)
                for i in range(numPoints):
                    for j in range(numTags):
                        if radiiKeys[j] == str(radii[i]):
                            tags[i] = tagValues[j]

            elif tags is not None:
                msg = "Using tag data for colour information in "
                msg += "BallPlot.setData()"
                debugMsg(msg)

                # check that the number of tags is correct
                if len(tags) != numPoints:
                    errorString = "The number of tags needs to be the"
                    errorString += "same as the number of points"
                    raise ValueError, errorString

                # need to find out the number of different tags
                valueDict = {}
                for i in range(numPoints):
                    valueDict[tags[i]] = 1
                numTags = len(valueDict.keys())
                tagValues = valueDict.keys()
                tagValues.sort()

            # now scale the tags
            scaledTags = Numeric.zeros(numPoints, typecode=Numeric.Float)
            if numTags == 1:
                pass
            else:
                for i in range(numPoints):
                    for j in range(numTags):
                        if tagValues[j] == tags[i]:
                            scaledTags[i] = float(j)/float(numTags-1)
        else:
            # barf
            raise ValueError, \
                    "Cannot construct BallPlot with the given input.  Exiting."

        # share the data around
        self.x = x
        self.y = y
        self.z = z
        self.radii = radii

        # share the colours around
        self.red = red
        self.green = green
        self.blue = blue

        self.scene.centre = centre
        self.scene.bounds = bounds

        return

    def render(self):
        """
        Does BallPlot specific rendering tasks
        """
        debugMsg("Called render() in BallPlot")
        self.renderer.runString("// BallPlot.render()")

        # grab the data
        x = self.x
        y = self.y
        z = self.z
        radii = self.radii

        # grab the colours
        red = self.red
        green = self.green
        blue = self.blue

        evalString = ""
        for i in range(len(x)):
            evalString += "sphere {\n"
            evalString += "  <%f, %f, %f> %f\n" % (x[i], y[i], z[i], radii[i])
            evalString += "  pigment {\n"
            evalString += "    rgb <%f, %f, %f>\n" % (red[i], green[i], blue[i])
            evalString += "  }\n"
            evalString += "}\n"

        self.renderer.runString(evalString)

        # set the title if set
        if self.title is not None:
            # not implemented yet
            pass

        return
 

class ContourPlot(Plot):
    """
    Contour plot
    """
    def __init__(self, scene):
        """
        Initialisation of the ContourPlot class
        
        @param scene: the scene with which to associate the contour plot
        @type scene: Scene object
        """
        debugMsg("Called ContourPlot.__init__()")
        Plot.__init__()

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to the plot

        @param dataList: list/tuple of data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in ContourPlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"
        
        return

class LinePlot(Plot):
    """
    Line plot
    """
    def __init__(self, scene):
        """
        Initialisation of the ContourPlot class
        
        @param scene: the scene with which to associate the line plot
        @type scene: Scene object
        """
        debugMsg("Called LinePlot.__init__()")
        Plot.__init__()

        if scene is None:
            raise ValueError, "You must specify a scene object"

    def setData(self, *dataList):
        """
        Set data to the plot

        @param dataList: list/tuple of data to set to the plot
        @type dataList: tuple
        """
        debugMsg("Called setData() in LinePlot()")

        if dataList is None:
            raise ValueError, "You must specify a data list"
        
        return

# vim: expandtab shiftwidth=4:

