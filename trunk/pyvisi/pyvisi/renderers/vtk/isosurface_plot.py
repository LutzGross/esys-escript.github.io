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

# $Id: isosurface_plot.py,v 1.1 2005/11/30 03:07:18 paultcochrane Exp $

"""
Class and functions associated with a pyvisi IsosurfacePlot objects
"""

# generic imports
from pyvisi.renderers.vtk.common import debugMsg
import Numeric
import os
import copy

# module specific imports
from pyvisi.renderers.vtk.plot import Plot

__revision__ = '$Revision: 1.1 $'

class IsosurfacePlot(Plot):
    """
    Isosurface plot
    """
    def __init__(self, scene):
        """
        Initialisation of the IsosurfacePlot class
        
        @param scene: The Scene to render the plot in
        @type scene: Scene object
        """
        debugMsg("Called IsosurfacePlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer
        self.renderer.addToInitStack("# IsosurfacePlot.__init__()")

        # labels and stuff
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None

        # how many contours?
        self.numContours = 5

        # contour range
        self.contMin = None
        self.contMax = None

        # default values for fname, format and scalars
        self.fname = None
        self.format = None
	self.scalars = None

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

	@param scalars: the name of the scalar data in the vtk file to use
	@type scalars: string
        """
        debugMsg("Called setData() in IsosurfacePlot()")

        self.renderer.runString("# IsosurfacePlot.setData()")

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
	## scalars
	if options.has_key('scalars'):
	    scalars = options['scalars']
	else:
	    scalars = None

	# we want to pass this info around
	self.fname = fname
	self.format = format
	self.scalars = scalars

	# reset the default values for shared info
	self.escriptData = False
	self.otherData = False

        # do a sanity check on the inputs
	if len(dataList) == 0 and fname is None:
	    raise ValueError, \
		    "You must specify a data list or an input filename"

	if len(dataList) != 0 and fname is not None:
	    raise ValueError, \
		    "You cannot specify a data list as well as an input file"

	if fname is not None and scalars is None:
	    debugMsg("No scalars specified; using default in vtk")

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

	    # now check the dimensionality of the data and the grid; make
	    # sure everything looks right before doing anything more

	    raise ImplementationError, "Can't process escript Data yet"

	elif self.otherData:

	    raise ImplementationError, "Can't process plain array data yet"

	# run the stuff for when we're reading from file
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

	    # need to do a delaunay 3D here to get decent looking isosurfaces
	    evalString = "_del3D = vtk.vtkDelaunay3D()\n"
	    evalString += "_del3D.SetInput(_reader.GetOutput())\n"
	    evalString += "_del3D.SetOffset(2.5)\n"
	    evalString += "_del3D.SetTolerance(0.001)\n"
	    evalString += "_del3D.SetAlpha(0.0)"
	    self.renderer.runString(evalString)

	    # get the model centre and bounds
	    evalString = "_centre = _reader.GetOutput().GetCenter()\n"
	    evalString += "_bounds = _reader.GetOutput().GetBounds()"
	    self.renderer.runString(evalString)

    def render(self):
        """
        Does IsosurfacePlot object specific (pre)rendering stuff
        """
        debugMsg("Called IsosurfacePlot.render()")

        self.renderer.runString("# IsosurfacePlot.render()")

        # set up a contour filter
        evalString = "_cont = vtk.vtkContourGrid()\n"
        evalString += "_cont.SetInput(_del3D.GetOutput())\n"

        # if contMin and contMax are or aren't set then handle the different
        # situations
        if self.contMin is not None and self.contMax is not None:
            evalString += "_cont.GenerateValues(%d, %f, %f)\n" % \
                    (self.numContours, self.contMin, self.contMax)
        elif self.contMin is not None and self.contMax is None:
            evalString += "(_contMin, _contMax) = _reader.GetOutput()."
            evalString += "GetPointData().GetScalars().GetRange()\n"
            evalString += "_cont.GenerateValues(%d, %f, _contMax)\n" % \
                    (self.numContours, self.contMin)
        elif self.contMin is None and self.contMax is not None:
            evalString += "(_contMin, _contMax) = _reader.GetOutput()."
            evalString += "GetPointData().GetScalars().GetRange()\n"
            evalString += "_cont.GenerateValues(%d, _contMin, %f)\n" % \
                    (self.numContours, self.contMax)
        elif self.contMin is None and self.contMax is None:
            evalString += "(_contMin, _contMax) = _reader.GetOutput()."
            evalString += "GetPointData().GetScalars().GetRange()\n"
            evalString += "_cont.GenerateValues(%d, _contMin, _contMax)\n" % \
                    (self.numContours)
        else:
            # barf, really shouldn't have got here
            raise ValueError, \
                    "Major problems in IsosurfacePlot: contMin and contMax"

        evalString += "_cont.GenerateValues(5, 0.25, 0.75)\n"
        evalString += "_cont.ComputeScalarsOn()"
        self.renderer.runString(evalString)

        # set up the mapper
        evalString = "_mapper = vtk.vtkDataSetMapper()\n"
        evalString += "_mapper.SetInput(_cont.GetOutput())\n"
        evalString += "_mapper.ScalarVisibilityOn()"
        self.renderer.runString(evalString)

        # set up the actor
        evalString = "_actor = vtk.vtkActor()\n"
        evalString += "_actor.SetMapper(_mapper)"
        self.renderer.runString(evalString)

        # add to the renderer
        evalString = "_renderer.AddActor(_actor)"
        self.renderer.runString(evalString)

        # set up the text properties for nice text
        evalString = "_font_size = 18\n"
        evalString += "_textProp = vtk.vtkTextProperty()\n"
        evalString += "_textProp.SetFontSize(_font_size)\n"
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

        # put an outline around the data
        evalString = "_outline = vtk.vtkOutlineSource()\n"
        evalString += "_outline.SetBounds(_bounds)\n"

        # make its mapper
        evalString += "_outlineMapper = vtk.vtkPolyDataMapper()\n"
        evalString += "_outlineMapper.SetInput(_outline.GetOutput())\n"

        # make its actor
        evalString += "_outlineActor = vtk.vtkActor()\n"
        evalString += "_outlineActor.SetMapper(_outlineMapper)\n"
        evalString += "_outlineActor.GetProperty().SetColor(0,0,0)"
        self.renderer.runString(evalString)

        # add to the renderer
        evalString = "_renderer.AddActor(_outlineActor)"
        self.renderer.runString(evalString)

        # make a lookup table for the colour map and invert it (colours look
        # better when it's inverted)
        evalString = "_lut = vtk.vtkLookupTable()\n"
        evalString += "_refLut = vtk.vtkLookupTable()\n"
        evalString += "_lut.Build()\n"
        evalString += "_refLut.Build()\n"
        evalString += "for _j in range(256):\n"
        evalString += "    _lut.SetTableValue(_j, "
        evalString += "_refLut.GetTableValue(255-_j))"
        self.renderer.runString(evalString)

        # add some axes
        evalString = "_axes = vtk.vtkCubeAxesActor2D()\n"
        evalString += "_axes.SetInput(_reader.GetOutput())\n"
        evalString += "_axes.SetCamera(_renderer.GetActiveCamera())\n"
        evalString += "_axes.SetLabelFormat(\"%6.4g\")\n"
        evalString += "_axes.SetFlyModeToOuterEdges()\n"
        evalString += "_axes.SetFontFactor(0.8)\n"
        evalString += "_axes.SetAxisTitleTextProperty(_textProp)\n"
        evalString += "_axes.SetAxisLabelTextProperty(_textProp)\n"
        ### xlabel
        if self.xlabel is not None:
            evalString += "_axes.SetXLabel(\"%s\")\n" % self.xlabel
        else:
            evalString += "_axes.SetXLabel(\"\")\n"
        ### ylabel
        if self.ylabel is not None:
            evalString += "_axes.SetYLabel(\"%s\")\n" % self.ylabel
        else:
            evalString += "_axes.SetYLabel(\"\")\n"
        ### zlabel
        if self.zlabel is not None:
            evalString += "_axes.SetZLabel(\"%s\")\n" % self.zlabel
        else:
            evalString += "_axes.SetZLabel(\"\")\n"
        evalString += "_axes.SetNumberOfLabels(5)\n"
        evalString += "_axes.GetProperty().SetColor(0,0,0)"
        self.renderer.runString(evalString)

        # add to the renderer
        evalString = "_renderer.AddActor(_axes)"
        self.renderer.runString(evalString)

        # play around with lighting
        evalString = "_light1 = vtk.vtkLight()\n"
        evalString += "_light1.SetFocalPoint(_centre)\n"
        evalString += "_light1.SetPosition(_centre[0]-_bounds[1], "
        evalString += "_centre[1]-_bounds[3], _centre[2]+_bounds[5])\n"
        evalString += "_renderer.AddLight(_light1)\n"
        evalString += "_light2 = vtk.vtkLight()\n"
        evalString += "_light2.SetFocalPoint(_centre)\n"
        evalString += "_light2.SetPosition(_centre[0]+_bounds[1], "
        evalString += "_centre[1]+_bounds[3], _centre[2]-_bounds[5])\n"
        evalString += "_renderer.AddLight(_light2)"
        self.renderer.runString(evalString)

        # this shouldn't be here!!!!
        self.renderer.runString("_renderer.SetBackground(1,1,1)")

        return


# vim: expandtab shiftwidth=4:
