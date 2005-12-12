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

# $Id: ellipsoid_plot.py,v 1.1 2005/11/30 03:07:18 paultcochrane Exp $

"""
Class and functions associated with a pyvisi EllipsoidPlot objects
"""

# generic imports
from pyvisi.renderers.vtk.common import debugMsg
import Numeric
import os
import copy

# module specific imports
from pyvisi.renderers.vtk.plot import Plot

__revision__ = '$Revision: 1.1 $'

class EllipsoidPlot(Plot):
    """
    Ellipsoid plot
    """
    def __init__(self, scene):
        """
        Initialisation of the EllipsoidPlot class

        @param scene: The Scene to render the plot in
        @type scene: Scene object
        """
        debugMsg("Called EllipsoidPlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer
        self.renderer.addToInitStack("# EllipsoidPlot.__init__()")

        # labels and stuff
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None
        
        # default values for fname, format and tensors
        self.fname = None
        self.format = None
	self.tensors = None

	# default values for shared info
	self.escriptData = False
	self.otherData = False

        # add the plot to the scene
        scene.add(self)

    def setData(self, *dataList, **options):
        """
        Set data to the plot

        @param dataList: List of data to set to the plot
        @type dataList: tuple

        @param options: Dictionary of keyword options to the method
        @type options: dict

	@param fname: the name of the input vtk file
	@type fname: string

	@param format: the format of the input vtk file ('vtk' or 'vtk-xml')
	@type format: string

	@param tensors: the name of the tensor data in the vtk file to use
	@type tensors: string
        """
        debugMsg("Called setData() in EllipsoidPlot()")

        self.renderer.runString("# EllipsoidPlot.setData()")

        # process the options, if any
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
	## tensors
	if options.has_key('tensors'):
	    tensors = options['tensors']
	else:
	    tensors = None

        # do a quick sanity check on the inputs
        if fname is None or format is None:
            raise ValueError, "You must supply an input file and its format"

        # we want to pass this info around
        self.fname = fname
        self.format = format
	self.tensors = tensors

	# reset the default values for shared info
	self.escriptData = False
	self.otherData = False

	# do some sanity checking on the input args
	if len(dataList) == 0 and fname is None:
	    raise ValueError, \
		    "You must specify a data list or an input filename"

	if len(dataList) != 0 and fname is not None:
	    raise ValueError, \
		    "You cannot specify a data list as well as an input file"

	if fname is not None and tensors is None:
	    debugMsg("No tensors specified; using default in vtk")

	if fname is not None and format is None:
	    raise ValueError, "You must specify an input file format"

	if fname is None and format is not None:
	    raise ValueError, "Format specified but no input filename"

	# if have just a data list, check the objects passed in to see if
	# they are escript data objects or not
	if len(dataList) != 0:
	    for obj in dataList:
		try:
		    obj.convertToNumArray()
		    # ok, we've got escript data, set the flag
		    self.escriptData = True
		except AttributeError:
		    self.otherData = True

	# if we have both escript data and other data, barf as can't handle
	# that yet
	if self.escriptData and self.otherData:
	    raise TypeError, \
		    "Sorry can't handle both escript and other data yet"

	# now generate the code for the case when we have just escript data
	# passed into setData()
	if self.escriptData:
	    # get the relevant bits of data
	    if len(dataList) == 1:
		# only one data variable, will need to get the domain from it
		escriptZ = dataList[0]
		escriptX = escriptZ.getDomain().getX()
	    else:
		errorString = \
			"Expecting only 1 element in data list.  I got %d" \
			% len(dataList)
		raise ValueError, errorString

	    ####!!!! now process the data properly so that I can plot it
	    print "escriptZ.shape() = %s" % escriptZ.shape()
	    print "escriptX.shape() = %s" % escriptX.shape()

	    # need to check the dimensionality of the data, i.e. is it 2x2
	    # or 3x3 or what?  This will determine the code that gets
	    # generated.

	    raise ImplementationError, "Can't process escript Data yet"

	elif self.otherData:

	    # do some checks to make sure have the right kind of shape for
	    # the data and then generate the code

	    raise ImplementationError, "Can't process plain array data yet"

	if fname is not None:

	    # check to see if the file exists
	    if not os.path.exists(fname):
		raise SystemError, "File %s not found" % fname

	    if format == 'vtk':
		# read old-style vtk files
		evalString = "_reader = vtk.vtkUnstructuredGridReader()\n"
	    elif format == 'vtk-xml':
		# read vtk xml files
		evalString = "_reader = vtk.vtkXMLUnstructuredGridReader()\n"
	    else:
		# barf
		raise ValueError, "Unknown format.  I got %s" % format

	    evalString += "_reader.SetFileName(\"%s\")\n" % fname
	    evalString += "_reader.Update()"

	    self.renderer.runString(evalString)

	    # grab the grid of the data
	    self.renderer.runString("_grid = _reader.GetOutput()")

	    # convert the cell data to point data
	    evalString = "_c2p = vtk.vtkCellDataToPointData()\n"
	    evalString += "_c2p.SetInput(_grid)"
	    self.renderer.runString(evalString)

	    # now extract the tensor components
	    evalString = "_extract = vtk.vtkExtractTensorComponents()\n"
	    evalString += "_extract.SetInput(_c2p.GetOutput())\n"
	    evalString += "_extract.SetScalarModeToEffectiveStress()\n"
	    evalString += "_extract.ExtractScalarsOn()\n"
	    evalString += "_extract.PassTensorsToOutputOn()\n"
	    evalString += "_extract.ScalarIsEffectiveStress()\n"

	    evalString += "_extractGrid = _extract.GetOutput()\n"
	    evalString += "_extractGrid.Update()\n"
	    evalString += "_extractScalarRange = "
	    evalString += \
		    "_extractGrid.GetPointData().GetScalars().GetRange()\n"
	    self.renderer.runString(evalString)

        return

    def render(self):
        """
        Does EllipsoidPlot object specific (pre)rendering stuff
        """
        debugMsg("Called EllipsoidPlot.render()")

        self.renderer.runString("# EllipsoidPlot.render()")

        # make a sphere source for the glyphs
        evalString = "_sphere = vtk.vtkSphereSource()\n"
        evalString += "_sphere.SetThetaResolution(6)\n"
        evalString += "_sphere.SetPhiResolution(6)\n"
        evalString += "_sphere.SetRadius(0.5)"
        self.renderer.runString(evalString)

        # make tensor glyphs
        evalString = "_glyph = vtk.vtkTensorGlyph()\n"
        evalString += "_glyph.SetSource(_sphere.GetOutput())\n"
        evalString += "_glyph.SetInput(_extractGrid)\n"
        evalString += "_glyph.SetColorModeToScalars()\n"
        evalString += "_glyph.ScalingOn()\n"
        evalString += "_glyph.SetMaxScaleFactor(5.0)\n"
        evalString += "_glyph.SetScaleFactor(1.0)\n"
        evalString += "_glyph.ClampScalingOn()"
        self.renderer.runString(evalString)

        # make a stripper for faster rendering
        evalString = "_stripper = vtk.vtkStripper()\n"
        evalString += "_stripper.SetInput(_glyph.GetOutput())"
        self.renderer.runString(evalString)

        # make the normals of the data
        evalString = "_normals = vtk.vtkPolyDataNormals()\n"
        evalString += "_normals.SetInput(_stripper.GetOutput())"
        self.renderer.runString(evalString)

        # make the mapper for the data
        evalString = "_mapper = vtk.vtkPolyDataMapper()\n"
        evalString += "_mapper.SetInput(_normals.GetOutput())\n"
        evalString += "_mapper.SetLookupTable(_lut)\n"
        evalString += "_mapper.SetScalarRange(_extractScalarRange)"
        self.renderer.runString(evalString)

        # make the actor
        evalString = "_actor = vtk.vtkActor()\n"
        evalString += "_actor.SetMapper(_mapper)"
        self.renderer.runString(evalString)

        # add the actor
        self.renderer.runString("_renderer.AddActor(_actor)")

        # set up the text properties for nice text
        evalString = "_textProp = vtk.vtkTextProperty()\n"
        evalString += "_textProp.SetFontFamilyToArial()\n"
        evalString += "_textProp.BoldOff()\n"
        evalString += "_textProp.ItalicOff()\n"
        evalString += "_textProp.ShadowOff()\n"
        evalString += "_textProp.SetColor(0.0, 0.0, 0.0)"
        self.renderer.runString(evalString)

        # if a title is set, put it in here
        if self.title is not None:
            # make a title
            evalString = "_title = vtk.vtkTextMapper()\n"
            evalString += "_title.SetInput(\"%s\")\n" % self.title

            # make the title text use the text properties
            evalString += "_titleProp = _title.GetTextProperty()\n"
            evalString += "_titleProp.ShallowCopy(_textProp)\n"
            evalString += "_titleProp.SetJustificationToCentered()\n"
            evalString += "_titleProp.SetVerticalJustificationToTop()\n"
            evalString += "_titleProp.SetFontSize(20)\n"
            evalString += "_titleProp.BoldOn()\n"

            # make the actor for the title
            evalString += "_titleActor = vtk.vtkTextActor()\n"
            evalString += "_titleActor.SetMapper(_title)\n"
            evalString += "_titleActor.GetPositionCoordinate()."
            evalString += "SetCoordinateSystemToNormalizedDisplay()\n"
            evalString += "_titleActor.GetPositionCoordinate()."
            evalString += "SetValue(0.5, 0.95)"
            self.renderer.runString(evalString)

            # add to the renderer
            evalString = "_renderer.AddActor(_titleActor)"
            self.renderer.runString(evalString)

        # add a scalar bar
        evalString = "_scalarBar = vtk.vtkScalarBarActor()\n"
        evalString += "_scalarBar.SetLookupTable(_lut)\n"
        evalString += "_scalarBar.SetWidth(0.1)\n"
        evalString += "_scalarBar.SetHeight(0.8)\n"
        evalString += "_scalarBar.SetPosition(0.9, 0.15)"
        self.renderer.runString(evalString)

        # set up the label text properties 
        evalString = "_scalarBarTextProp = _scalarBar.GetLabelTextProperty()\n"
        evalString += "_scalarBarTextProp.ShallowCopy(_textProp)\n"
        evalString += "_scalarBarTextProp.SetFontSize(10)\n"

        evalString += "_renderer.AddActor(_scalarBar)"
        self.renderer.runString(evalString)

        return


# vim: expandtab shiftwidth=4:
