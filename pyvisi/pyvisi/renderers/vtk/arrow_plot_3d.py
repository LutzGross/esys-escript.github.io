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

# $Id: arrow_plot_3d.py,v 1.1 2005/11/30 03:07:18 paultcochrane Exp $

"""
Class and functions associated with a pyvisi ArrowPlot3D objects
"""

# generic imports
from pyvisi.renderers.vtk.common import debugMsg
import Numeric
import os
import copy

# module specific imports
from pyvisi.renderers.vtk.plot import Plot

__revision__ = '$Revision: 1.1 $'

class ArrowPlot3D(Plot):
    """
    Arrow field plot in three dimensions
    """
    def __init__(self, scene):
        """
        Initialisation of the ArrowPlot3D class
        
        @param scene: The Scene to render the plot in
        @type scene: Scene object
        """
        debugMsg("Called ArrowPlot3D.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer

        self.renderer.runString("# ArrowPlot3D.__init__()")

        # default values for fname, format and vectors
        self.fname = None
        self.format = None
	self.vectors = None
	
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

        @param options: Dictionary of extra options
        @type options: dict

        @param fname: Filename of the input vtk file
        @type fname: string

        @param format: Format of the input vtk file ('vtk' or 'vtk-xml')
        @type format: string

	@param vectors: the name of the vector data in the vtk file to use
	@type vectors: string
        """
        debugMsg("Called setData() in ArrowPlot3D()")
        self.renderer.runString("# ArrowPlot3D.setData()")

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
	## vectors
	if options.has_key('vectors'):
	    vectors = options['vectors']
	else:
	    vectors = None

        # we want to pass this info around
        self.fname = fname
        self.format = format
	self.vectors = vectors

	# reset the default values for shared info
	self.escriptData = False
	self.otherData = False

        # do some sanity checking on the inputs
	if len(dataList) == 0 and fname is None:
	    raise ValueError, \
		    "You must specify a data list or an input filename"

	if len(dataList) != 0 and fname is not None:
	    raise ValueError, \
		    "You cannot specify a data list as well as an input file"

	if fname is not None and vectors is None:
	    debugMsg("No vectors specified; using default in vtk")

	if fname is not None and format is None:
	    raise ValueError, "You must specify an input file format"

        if fname is None and format is not None:
            raise ValueError, "Format specified, but no input filename"

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

	# if we have both escript and other data, barf as can't handle that
	# just yet
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

	    # need to check the length of the vector data and the shape of
	    # the mesh to make sure is appropriate for an ArrowPlot3D
	    # object.

	    raise ImplementationError, "Can't process escript Data yet"

	elif self.otherData:
            # do some sanity checking on the data
            if len(dataList) != 6:
                errorStr = "Must have six vectors as input: x, y, z, "
                errorStr += "dx, dy, dz, found: %d" % len(dataList)
                raise ValueError, errorStr
                        
    
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
                xData = dataList[0]
                yData = dataList[1]
                zData = dataList[2]
                dxData = dataList[3]
                dyData = dataList[4]
                dzData = dataList[5]
            elif len(dataList[0].shape) == 2:
                xData = dataList[0].flat
                yData = dataList[1].flat
                zData = dataList[2].flat
                dxData = dataList[3].flat
                dyData = dataList[4].flat
                dzData = dataList[5].flat
            else:
                raise ValueError, "Input vectors can only be 1D or 2D"
    
            # now pass the data to the render dictionary so that the render
            # code knows what it's supposed to plot

            # x data
            self.renderer.renderDict['_x'] = copy.deepcopy(xData)
    
            # y data
            self.renderer.renderDict['_y'] = copy.deepcopy(yData)
    
            # z data
            self.renderer.renderDict['_z'] = copy.deepcopy(zData)
    
            # dx data
            self.renderer.renderDict['_dx'] = copy.deepcopy(dxData)
    
            # dy data
            self.renderer.renderDict['_dy'] = copy.deepcopy(dyData)
    
            # dz data
            self.renderer.renderDict['_dz'] = copy.deepcopy(dzData)
    
            # keep the number of points for future reference
            numPoints = len(xData)
    
            # construct the points data
            evalString = "_points = vtk.vtkPoints()\n"
            evalString += "_points.SetNumberOfPoints(%d)\n" % numPoints
            evalString += "for _j in range(%d):\n" % numPoints
            evalString += \
                    "    _points.InsertPoint(_j, _x[_j], _y[_j], _z[_j])\n"
            self.renderer.runString(evalString)
    
            # construct the vectors
            evalString = "_vectors = vtk.vtkFloatArray()\n"
            evalString += "_vectors.SetNumberOfComponents(3)\n"
            evalString += "_vectors.SetNumberOfTuples(%d)\n" % numPoints
            evalString += "_vectors.SetName(\"vectors\")\n"
            evalString += "for _j in range(%d):\n" % numPoints
            evalString += \
                    "    _vectors.InsertTuple3(_j, _dx[_j], _dy[_j], _dz[_j])\n"
            self.renderer.runString(evalString)
    
            # construct the grid
            evalString = "_grid = vtk.vtkUnstructuredGrid()\n"
            evalString += "_grid.SetPoints(_points)\n"
            evalString += "_grid.GetPointData().AddArray(_vectors)\n"
            evalString += "_grid.GetPointData().SetActiveVectors(\"vectors\")"
            self.renderer.runString(evalString)

	# run the stuff for when we're reading from file
	if fname is not None:

            # had best make sure it exists
            if not os.path.exists(fname):
                raise SystemError, "File: %s not found" % fname
            
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

            # grab the grid
            evalString = "_grid = _reader.GetOutput()"
            self.renderer.runString(evalString)

            # get the norm of the vectors
            evalString = "_norm = vtk.vtkVectorNorm()\n"
            evalString += "_norm.SetInput(_grid)"
            self.renderer.runString(evalString)

    def render(self):
        """
        Does ArrowPlot3D specific rendering tasks
        """
        debugMsg("Called render() in ArrowPlot3D")
        self.renderer.runString("# ArrowPlot3D.render()")

        # make the arrow source
        self.renderer.runString("_arrow = vtk.vtkArrowSource()")

        # get the maximum norm of the data
        evalString = "_maxNorm = _grid.GetPointData().GetVectors().GetMaxNorm()"
        self.renderer.runString(evalString)

        # make the glyph
        evalString = "_glyph = vtk.vtkGlyph3D()\n"
        evalString += "_glyph.ScalingOn()\n"
        evalString += "_glyph.SetScaleModeToScaleByVector()\n"
        evalString += "_glyph.SetColorModeToColorByVector()\n"
        evalString += "_glyph.SetScaleFactor(0.1/_maxNorm)\n"
        # if we have a vtk grid from file, use the norm output of that
        if self.fname is not None and self.format is not None:
            evalString += "_glyph.SetInput(_norm.GetOutput())\n"
        else:
            evalString += "_glyph.SetInput(_grid)\n"
        evalString += "_glyph.SetSource(_arrow.GetOutput())\n"
        evalString += "_glyph.ClampingOff()"
        self.renderer.runString(evalString)

        # set up a stripper for faster rendering
        evalString = "_stripper = vtk.vtkStripper()\n"
        evalString += "_stripper.SetInput(_glyph.GetOutput())"
        self.renderer.runString(evalString)

        # invert the lookup table, so the colours are nicer
        evalString = "_lut = vtk.vtkLookupTable()\n"
        evalString += "_lut.Build()\n"
        evalString += "_refLut = vtk.vtkLookupTable()\n"
        evalString += "_refLut.Build()\n"
        evalString += "for _i in range(256):\n"
        evalString += \
                "    _lut.SetTableValue(_i, _refLut.GetTableValue(255-_i))"
        self.renderer.runString(evalString)

        # set up the mapper
        evalString = "_mapper = vtk.vtkPolyDataMapper()\n"
        evalString += "_mapper.SetInput(_stripper.GetOutput())\n"
        evalString += "_mapper.SetScalarRange(0, _maxNorm)"
        self.renderer.runString(evalString)

        # set up the actor
        evalString = "_actor = vtk.vtkActor()\n"
        evalString += "_actor.SetMapper(_mapper)"
        self.renderer.runString(evalString)

        # add the actor
        self.renderer.runString("_renderer.AddActor(_actor)")

        # text properties
        evalString = "_font_size = 14\n"  # this will need to be an option!!
        evalString += "_textProp = vtk.vtkTextProperty()\n"
        evalString += "_textProp.SetFontSize(_font_size)\n"
        evalString += "_textProp.SetFontFamilyToArial()\n"
        evalString += "_textProp.BoldOff()\n"
        evalString += "_textProp.ItalicOff()\n"
        evalString += "_textProp.ShadowOff()\n"
        evalString += "_textProp.SetColor(0,0,0)\n"
        self.renderer.runString(evalString)

        # set the title if set
        if self.title is not None:
            # add a title
            evalString = "_titleMapper = vtk.vtkTextMapper()\n"
            evalString += "_titleMapper.SetInput(\"%s\")\n" % self.title
            
            evalString += "_titleProp = _titleMapper.GetTextProperty()\n"
            evalString += "_titleProp.ShallowCopy(_textProp)\n"
            evalString += "_titleProp.SetJustificationToCentered()\n"
            evalString += "_titleProp.SetVerticalJustificationToTop()\n"
            evalString += "_titleProp.SetFontSize(18)\n"
            
            # set up the text actor
            evalString += "_titleActor = vtk.vtkTextActor()\n"
            evalString += "_titleActor.SetMapper(_titleMapper)\n"
            evalString += "_titleActor.GetPositionCoordinate()."
            evalString += "SetCoordinateSystemToNormalizedDisplay()\n"
            evalString += "_titleActor.GetPositionCoordinate()."
            evalString += "SetValue(0.5, 0.95)\n"

            evalString += "_renderer.AddActor(_titleActor)"
            self.renderer.runString(evalString)

        # set up some axes
        evalString = "_axes = vtk.vtkCubeAxesActor2D()\n"
        evalString += "_axes.SetCamera(_renderer.GetActiveCamera())\n"
        evalString += "_axes.SetFlyModeToOuterEdges()\n"
        if self.fname is not None and self.format is not None:
            evalString += "_axes.SetInput(_grid)\n"
        else:
            evalString += "_axes.SetBounds(min(_x)-_maxNorm, max(_x)+_maxNorm, "
            evalString += "min(_y)-_maxNorm, max(_y)+_maxNorm, "
            evalString += "min(_z)-_maxNorm, max(_z)+_maxNorm)\n"

        if self.xlabel is None:
            evalString += "_axes.SetXLabel(\"\")\n"
        else:
            evalString += "_axes.SetXLabel(\"%s\")\n" % self.xlabel

        if self.ylabel is None:
            evalString += "_axes.SetYLabel(\"\")\n"
        else:
            evalString += "_axes.SetYLabel(\"%s\")\n" % self.ylabel

        if self.zlabel is None:
            evalString += "_axes.SetZLabel(\"\")\n"
        else:
            evalString += "_axes.SetZLabel(\"%s\")\n" % self.zlabel

        # set up the axes properties
        evalString += "_axesProp = _axes.GetProperty()\n"
        evalString += "_axesProp.SetColor(0,0,0)\n"

        # set up the axes title properties
        evalString += "_axesTitleProp = _axes.GetAxisTitleTextProperty()\n"
        evalString += "_axesTitleProp.ShallowCopy(_textProp)\n"
        
        # set up the axes label properties
        evalString += "_axesLabelProp = _axes.GetAxisLabelTextProperty()\n"
        evalString += "_axesLabelProp.ShallowCopy(_textProp)\n"
        evalString += "_axesLabelProp.SetFontSize(8)\n"
        self.renderer.runString(evalString)

        # add the axes to the renderer
        self.renderer.runString("_renderer.AddActor(_axes)")

        # reset the camera, will make things look nicer
        ### is this the right place to put this???
        self.renderer.runString("_renderer.ResetCamera()")
        
        ### this should be somewhere else too...
        self.renderer.runString("_renderer.SetBackground(1,1,1)")

        return

# vim: expandtab shiftwidth=4:
