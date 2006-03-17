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

# $Id: ball_plot.py,v 1.1 2005/11/30 03:07:18 paultcochrane Exp $

"""
Class and functions associated with a pyvisi BallPlot objects
"""

# generic imports
from pyvisi.renderers.vtk.common import debugMsg
import Numeric
import os
import copy

# module specific imports
from pyvisi.renderers.vtk.plot import Plot

__revision__ = '$Revision: 1.1 $'

class BallPlot(Plot):
    """
    Ball plot
    """
    def __init__(self, scene):
        debugMsg("Called BallPlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer

        self.renderer.runString("# BallPlot.__init__()")

        # add the plot to the scene
        scene.add(self)

    def setData(self, points=None, 
            fname=None, format=None,
            radii=None, colors=None, tags=None):
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
        self.renderer.runString("# BallPlot.setData()")

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
            if radii is None:
                raise ValueError, \
                "You must specify the name of the scalars to use as the radius"
            if colors is None and tags is None:
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
	    # check to make sure the file exists
	    if not os.path.exists(fname):
		raise SystemError, "File: %s not found" % fname

            if format == "vtk-xml":
                debugMsg("Using vtk-xml file as input")
                # create the reader of the file
                evalString = "_reader = vtk.vtkXMLUnstructuredGridReader()\n"
            elif format == "vtk":
                debugMsg("Using old-style vtk file as input")
                # create the reader of the file
                evalString = "_reader = vtk.vtkUnstructuredGridReader()\n"

            evalString += "_reader.SetFileName(\"%s\")\n" % fname
            evalString += "_reader.Update()"
            self.renderer.runString(evalString)

            # read the output to an unstructured grid
            self.renderer.runString("_grid = _reader.GetOutput()")

            # note that these next few steps are only necessary in vtk 4.2,
            # 4.4 grab the data to use for the radii of the balls
            evalString = \
                    "_radii = _grid.GetPointData().GetScalars(\"%s\")" % \
                    radii
            self.renderer.runString(evalString)
    
            # grab the data to use for colouring the balls
            if colors is None:
                evalString = \
                    "_colours = _grid.GetPointData().GetScalars(\"%s\")" %\
                    tags
            else:
                evalString = \
                    "_colours = _grid.GetPointData().GetScalars(\"%s\")" % \
                    colors
            self.renderer.runString(evalString)
    
            # now work out the number of tags, and their values
            evalString = "_numPoints = _colours.GetNumberOfTuples()\n"
            evalString += "_valueDict = {}\n"
            evalString += "for _i in range(_numPoints):\n"
            evalString += "    _colourValue = _colours.GetValue(_i)\n"
            evalString += "    _valueDict[_colourValue] = 1\n"
    
            evalString += "_numColours = len(_valueDict.keys())\n"
    
            evalString += "_colourValues = _valueDict.keys()\n"
            evalString += "_colourValues.sort()"
            self.renderer.runString(evalString)
    
            # now count the number of colours, and make an evenly spaced
            # array of points between zero and one, then use these as the
            # scalars to colour by
            evalString = "_scaledColours = vtk.vtkFloatArray()\n"
            evalString += "_scaledColours.SetNumberOfTuples(_numPoints)\n"
            evalString += "_scaledColours.SetNumberOfComponents(1)\n"
            evalString += "_scaledColours.SetName(\"scaledColours\")\n"
            evalString += "for _i in range(_numPoints):\n"
            evalString += "    _colourValue = _colours.GetValue(_i)\n"
            evalString += "    for _j in range(_numColours):\n"
            evalString += "        if _colourValues[_j] == _colourValue:\n"
            evalString += "            _scaledColours.InsertTuple1(_i,"
            evalString += "float(_j)/float(_numColours-1))"
            self.renderer.runString(evalString)
    
            # now set up an array of two components to get the data through
            # the glyph object to the mapper (this is so that colouring and
            # scalaing work properly)
            evalString = "_data = vtk.vtkFloatArray()\n"
            evalString += "_data.SetNumberOfComponents(3)\n"
            evalString += \
                    "_data.SetNumberOfTuples(_radii.GetNumberOfTuples())\n"
            evalString += "_data.CopyComponent(0, _radii, 0)\n"
            evalString += "_data.CopyComponent(1, _colours, 0)\n"
            evalString += "_data.CopyComponent(2, _scaledColours, 0)\n"
            evalString += "_data.SetName(\"data\")"
            self.renderer.runString(evalString)
    
            # add the data array to the grid
            evalString = "_grid.GetPointData().AddArray(_data)\n"
    
            # make the data the active scalars
            evalString += "_grid.GetPointData().SetActiveScalars(\"data\")"
            self.renderer.runString(evalString)

        elif format is None and points is not None:
            debugMsg("Using user-defined point data in BallPlot.setData()")
            ### if we get to here, then we have to construct the points,
            ### the radii and the colours all from scratch, add them to the
            ### grid and get everything in the same form that we have for the
            ### case where we load a vtk data file

            numPoints = len(points)
            numRadii = len(radii)
            # do some sanity checking on the data
            if numRadii != numPoints:
                raise ValueError, \
                    "The number of points does not equal the number of radii"

            ### construct the grid from the point data
            self.renderer.renderDict['_pointData'] = copy.deepcopy(points)
            self.renderer.renderDict['_radiiData'] = copy.deepcopy(radii)
            # make the points
            evalString = "_points = vtk.vtkPoints()\n"
            evalString += "_points.SetNumberOfPoints(%d)\n" % numPoints
            evalString += "for _j in range(%d):\n" % numPoints
            evalString += "    _point = _pointData[_j]\n"
            evalString += \
                "    _points.InsertPoint(_j,_point[0],_point[1],_point[2])\n"
            self.renderer.runString(evalString)

            # make the radii
            evalString = "_radii = vtk.vtkFloatArray()\n"
            evalString += "_radii.SetNumberOfComponents(1)\n"
            evalString += "_radii.SetNumberOfValues(%d)\n" % numPoints
            evalString += "for _j in range(%d):\n" % numPoints
            evalString += "    _radii.InsertValue(_j, _radiiData[_j])\n"
            self.renderer.runString(evalString)

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

            self.renderer.renderDict['_tagData'] = copy.deepcopy(tags)

            # give the tag data to vtk
            evalString = "_tags = vtk.vtkFloatArray()\n"
            evalString += "_tags.SetNumberOfValues(%d)\n" % numPoints
            evalString += "_tags.SetNumberOfComponents(1)\n"
            evalString += "_tags.SetName(\"tags\")\n"
            evalString += "for _j in range(%d):\n" % numPoints
            evalString += "    _tags.InsertValue(_j, _tagData[_j])\n"
            self.renderer.runString(evalString)

            # now scale the tags
            scaledTags = Numeric.zeros(numPoints, typecode=Numeric.Float)
            if numTags == 1:
                pass
            else:
                for i in range(numPoints):
                    for j in range(numTags):
                        if tagValues[j] == tags[i]:
                            scaledTags[i] = float(j)/float(numTags-1)

            self.renderer.renderDict['_scaledTagData'] = \
		    copy.deepcopy(scaledTags)

            # now give vtk the scaled tag data
            evalString = "_scaledTags = vtk.vtkFloatArray()\n"
            evalString += "_scaledTags.SetNumberOfValues(%d)\n" % numPoints
            evalString += "_scaledTags.SetNumberOfComponents(1)\n"
            evalString += "_scaledTags.SetName(\"scaledTags\")\n"
            evalString += "for _j in range(%d):\n" % numPoints
            evalString += \
                    "    _scaledTags.InsertValue(_j, _scaledTagData[_j])\n"
            self.renderer.runString(evalString)

            # now construct the data array
            ### this is a vtk 4.2, 4.4 specific thing.  vtk 4.5 and above
            ### have a better way to do it, but this is here for backwards 
            ### compatibility
            evalString = "_data = vtk.vtkFloatArray()\n"
            evalString += "_data.SetNumberOfComponents(3)\n"
            evalString += \
                    "_data.SetNumberOfTuples(_radii.GetNumberOfTuples())\n"
            evalString += "_data.CopyComponent(0, _radii, 0)\n"
            evalString += "_data.CopyComponent(1, _tags, 0)\n"
            evalString += "_data.CopyComponent(2, _scaledTags, 0)\n"
            evalString += "_data.SetName(\"data\")\n"
            self.renderer.runString(evalString)

            # now construct the grid
            evalString = "_grid = vtk.vtkUnstructuredGrid()\n"
            evalString += "_grid.SetPoints(_points)\n"

            # add the data array to the grid
            evalString += "_grid.GetPointData().AddArray(_data)\n"

            # make the data the active scalars
            evalString += "_grid.GetPointData().SetActiveScalars(\"data\")\n"

            self.renderer.runString(evalString)
        else:
            # barf
            raise ValueError, \
                    "Cannot construct BallPlot with the given input.  Exiting."

        return

    def render(self):
        """
        Does BallPlot specific rendering tasks
        """
        debugMsg("Called render() in BallPlot")
        self.renderer.runString("# BallPlot.render()")

        # to make sphere glyphs need a sphere source
        evalString = "_sphere = vtk.vtkSphereSource()\n"
        evalString += "_sphere.SetRadius(1.0)\n"
        evalString += "_sphere.SetThetaResolution(5)\n"
        evalString += "_sphere.SetPhiResolution(5)"
        self.renderer.runString(evalString)

        # the spheres are 3D glyphs so set that up
        evalString = "_glyph = vtk.vtkGlyph3D()\n"
        evalString += "_glyph.ScalingOn()\n"
        evalString += "_glyph.SetScaleModeToScaleByScalar()\n"
        evalString += "_glyph.SetColorModeToColorByScalar()\n"
        evalString += "_glyph.SetScaleFactor(1.0)\n"
        evalString += "_glyph.SetInput(_grid)\n"
        evalString += "_glyph.SetSource(_sphere.GetOutput())\n"
        evalString += "_glyph.ClampingOff()"
        self.renderer.runString(evalString)

        # set up a stripper (this will speed up rendering)
        evalString = "_stripper = vtk.vtkStripper()\n"
        evalString += "_stripper.SetInput(_glyph.GetOutput())\n"

        # denote the stripper as being before the mapper by default, and let
        # subsequent objects redefine this if necessary
        evalString += "_preMapper = _stripper"
        self.renderer.runString(evalString)

        # if any clip objects etc are registered, then get them to render
        # themselves here
        for obj in self.objectList:
            obj.render()

        # set up the mapper
        evalString = "_mapper = vtk.vtkPolyDataMapper()\n"
        evalString += "_mapper.SetInput(_preMapper.GetOutput())\n"
        evalString += "_mapper.ScalarVisibilityOn()\n"
        # note: this is for vtk 4.2, 4.4 (4.5 and above have a better
        # technique to colour the scalars, but that version isn't yet
        # standard, or in fact released)
        evalString += "_mapper.ColorByArrayComponent(\"data\", 2)\n"
        # should be done in setData()
        evalString += "_mapper.SetScalarRange(0, 1)"
        self.renderer.runString(evalString)

        # set up the actor
        evalString = "_actor = vtk.vtkActor()\n"
        evalString += "_actor.SetMapper(_mapper)"
        self.renderer.runString(evalString)

        # add the actor to the scene
        self.renderer.runString("_renderer.AddActor(_actor)")

        # set the title if set
        if self.title is not None:
            # text properties
            evalString = "_font_size = 20\n"  # this will need to be an option!!
            evalString += "_textProp = vtk.vtkTextProperty()\n"
            evalString += "_textProp.SetFontSize(_font_size)\n"
            evalString += "_textProp.SetFontFamilyToArial()\n"
            evalString += "_textProp.BoldOff()\n"
            evalString += "_textProp.ItalicOff()\n"
            evalString += "_textProp.ShadowOff()\n"
            evalString += "_textProp.SetColor(0,0,0)\n"
        
            # add a title
            evalString += "_titleMapper = vtk.vtkTextMapper()\n"
            evalString += "_titleMapper.SetInput(\"%s\")\n" % self.title
            
            evalString += "_titleProp = _titleMapper.GetTextProperty()\n"
            evalString += "_titleProp.ShallowCopy(_textProp)\n"
            evalString += "_titleProp.SetJustificationToCentered()\n"
            evalString += "_titleProp.SetVerticalJustificationToTop()\n"
            evalString += "_titleProp.BoldOn()\n"
            
            # set up the text actor
            evalString += "_titleActor = vtk.vtkTextActor()\n"
            evalString += "_titleActor.SetMapper(_titleMapper)\n"
            evalString += "_titleActor.GetPositionCoordinate()."
            evalString += "SetCoordinateSystemToNormalizedDisplay()\n"
            evalString += "_titleActor.GetPositionCoordinate()."
            evalString += "SetValue(0.5, 0.95)\n"

            evalString += "_renderer.AddActor(_titleActor)"
            self.renderer.runString(evalString)

        return
 
# vim: expandtab shiftwidth=4:
