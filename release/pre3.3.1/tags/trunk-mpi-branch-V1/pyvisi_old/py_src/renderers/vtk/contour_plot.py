"""
Class and functions associated with a pyvisi ContourPlot objects

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


# generic imports
from common import debugMsg
import numarray
import os
import copy

# module specific imports
from plot import Plot

class ContourPlot(Plot):
    """
    Contour plot
    """
    def __init__(self, scene):
        """
        Initialisation of the ContourPlot class
        
        @param scene: The Scene to render the plot in
        @type scene: Scene object
        """
        debugMsg("Called ContourPlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer

        # default values for shared info
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

        @param options: Dictionary of extra options
        @type options: dict

        @keyword fname: the name of the input vtk file
        @type fname: string

        @keyword format: the format of the input vtk file ('vtk' or 'vtk-xml')
        @type format: string

        @keyword scalars: the scalar data in the vtk file to use
        @type scalars: string
        """
        debugMsg("Called setData() in ContourPlot()")
        self.renderer.runString("# ContourPlot.setData()")

        # get the options, if any
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

        # do some sanity checking on the input args
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
                    "Sorry, can't handle both escript and other data yet"
        
        elif self.escriptData and not self.otherData:
            # do we have access to escript??
            try:
                # escript objects should be able to be converted to numarrays
                # so use this as our test if escript is available
                dataList[0].convertToNumArray()
                debugMsg("Using escript")
            except AttributeError:
                raise ImportError, "Unable to use escript"

        else:
            # well, just try and handle the data normally
            pass

        # now generate the code for the case when we have escript data
        # just passed into setData()
        if self.escriptData:
            # now need to check if have Finley or Bruce mesh
            # note that Bruce isn't actually implemented yet
            # do I need to worry about this in vtk???  I'll be just using
            # an unstructured grid anyway...
            
            ##!!!!  need to check for rank of data so know if it is 
            ##!!!!  scalar or not, if other than scalar have to do
            ##!!!!  something other than am doing here

            # get the relevant bits of data
            if len(dataList) == 1:
                # only one data variable, will need to get the domain from it
                ### the capital letter denotes this is an escript object
                escriptZ = dataList[0]
                escriptX = escriptZ.getDomain().getX()
            elif len(dataList) == 2:
                # first variable should be the domain, the second the data
                escriptX = dataList[0]
                escriptZ = dataList[1]
            else:
                errorString = \
                        "Expecting 1 or 2 elements in data list.  I got: %d" \
                        % len(dataList)
                raise ValueError, errorString

            # convert the data to numarray
            xData = escriptX[0].convertToNumArray()
            yData = escriptX[1].convertToNumArray()
            zData = escriptZ.convertToNumArray()

            # pass the data through to the pyvisi renderer
            ### the x data
            self.renderer.renderDict['_x'] = copy.deepcopy(xData)

            ### the y data
            self.renderer.renderDict['_y'] = copy.deepcopy(yData)

            ### the z data
            self.renderer.renderDict['_z'] = copy.deepcopy(zData)

            # calculate the max and min of the z data
            evalString = "_zMin = min(_z)\n"
            evalString += "_zMax = max(_z)"
            self.renderer.runString(evalString)

            # create the points
            evalString = "_points = vtk.vtkPoints()\n"
            evalString += "_points.SetNumberOfPoints(len(_x))\n"
            evalString += "for _i in range(len(_x)):\n"
            evalString += "    _points.InsertPoint(_i, _x[_i], _y[_i], 0)"
            self.renderer.runString(evalString)

            # create the data
            evalString = "_data = vtk.vtkFloatArray()\n"
            evalString += "_data.SetNumberOfComponents(1)\n"
            evalString += "_data.SetNumberOfValues(len(_z))\n"
            evalString += "for _i in range(len(_z)):\n"
            evalString += "    _data.InsertValue(_i, _z[_i])"
            self.renderer.runString(evalString)

            # set up the grid (it's polydata since we're doing a Delaunay2D)
            evalString = "_grid = vtk.vtkPolyData()\n"
            evalString += "_grid.SetPoints(_points)\n"
            evalString += "_grid.GetPointData().SetScalars(_data)"
            self.renderer.runString(evalString)

        elif self.otherData:
            # in this case, we can only accept data lists of length 1 or 3
            # a length of 2 creates an abiguity
            if len(dataList) == 1:
                # need to autogenerate the x and y data
                zData = dataList[0]
                # check that the zData has the right shape
                if len(zData.shape) != 2:
                    raise ValueError, \
                            "z data array is not of correct shape: %s" % \
                            zData.shape
                # autogen the x and y data
                zShape = zData.shape
                xData = range(1, zShape[0]+1)
                yData = range(1, zShape[1]+1)
                # check was created correctly just in case
                if len(xData) != zShape[0]:
                    raise ValueError, \
                        "Autogenerated xData not equal to first dim of zData"
                if len(yData) != zShape[1]:
                    raise ValueError, \
                        "Autogenerated yData not equal to second dim of zData"

            elif len(dataList) == 2:
                raise ValueError, \
                    "The data list can't be of length 2 for non-escript data"
            elif len(dataList) == 3:
                # now just pass the x, y and z data through to the renderer
                xData = dataList[0]
                yData = dataList[1]
                zData = dataList[2]
            else:
                raise ValueError, \
                    "Expecting a data list length of 1 or 3.  I got: %d" \
                    % len(dataList)

            # check the shapes of the data
            if len(xData.shape) != 1:
                raise ValueError, "x data array is not of correct shape: %s"% \
                        xData.shape

            if len(yData.shape) != 1:
                raise ValueError, "y data array is not of correct shape: %s"% \
                        yData.shape

            if len(zData.shape) != 2:
                raise ValueError, "z data array is not of correct shape: %s"% \
                        zData.shape

            # stringify the data to then pass to the renderer
            ### x data
            self.renderer.renderDict['_x'] = copy.deepcopy(xData)

            ### y data
            self.renderer.renderDict['_y'] = copy.deepcopy(yData)

            ### z data
            self.renderer.renderDict['_z'] = copy.deepcopy(zData)

            # calculate the min and max
            evalString = "_zMax = max(_z.flat)\n"
            evalString += "_zMin = min(_z.flat)"
            self.renderer.runString(evalString)

            # create the points
            evalString = "_points = vtk.vtkPoints()\n"
            evalString += "_points.SetNumberOfPoints(len(_x)*len(_y))\n"
            evalString += "_count = 0\n"
            evalString += "for _i in range(len(_x)):\n"
            evalString += "  for _j in range(len(_y)):\n"
            evalString += "    _points.InsertPoint(_count, _x[_i], _y[_j], 0)\n"
            evalString += "    _count += 1"
            self.renderer.runString(evalString)

            # create the data
            evalString = "_data = vtk.vtkFloatArray()\n"
            evalString += "_data.SetNumberOfComponents(1)\n"
            evalString += "_data.SetNumberOfValues(len(_x)*len(_y))\n"
            evalString += "_count = 0\n"
            evalString += "for _i in range(len(_x)):\n"
            evalString += "  for _j in range(len(_y)):\n"
            evalString += "    _data.InsertValue(_count, _z[_i][_j])\n"
            evalString += "    _count += 1"
            self.renderer.runString(evalString)

            # set up the grid (it's polydata since we're doing a Delaunay2D)
            evalString = "_grid = vtk.vtkPolyData()\n"
            evalString += "_grid.SetPoints(_points)\n"
            evalString += "_grid.GetPointData().SetScalars(_data)"
            self.renderer.runString(evalString)

        # run the stuff for when we're reading from file
        if fname is not None:
	    # make sure the file exists
	    if not os.path.exists(fname):
		raise SystemError, "File: %s not found" % fname

	    # now handle the different kinds of vtk files
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
    
            # read the output input an unstructured grid
            evalString = "_grid = _reader.GetOutput()\n"
            evalString += \
                    "_grid.GetPointData().SetActiveScalars(\"%s\")" % scalars
            self.renderer.runString(evalString)
    
            # grab the range of scalars for appropriate scaling of the colourmap
            evalString = \
                "_scalarRange = _grid.GetPointData().GetScalars().GetRange()\n"
            evalString += "_scalarMin = _scalarRange[0]\n"
            evalString += "_scalarMax = _scalarRange[1]\n"
            self.renderer.runString(evalString)

        return

    def render(self):
        """
        Does ContourPlot object specific (pre)rendering stuff
        """
        debugMsg("Called ContourPlot.render()")

        self.renderer.runString("# ContourPlot.render()")

        # set up the lookup table and reverse the order of the colours
        evalString = "_lut = vtk.vtkLookupTable()\n"
        evalString += "_lut.Build()\n"
        evalString += "_refLut = vtk.vtkLookupTable()\n"
        evalString += "_refLut.Build()\n"
        evalString += "for _i in range(256):\n"
        evalString += \
                "    _lut.SetTableValue(_i, _refLut.GetTableValue(255-_i))"
        self.renderer.runString(evalString)

        if self.escriptData or self.otherData:
            # triangulate the data
            evalString = "_delaunay = vtk.vtkDelaunay2D()\n"
            evalString += "_delaunay.SetInput(_grid)\n"
            evalString += "_delaunay.SetTolerance(0.001)"
            self.renderer.runString(evalString)

            # set up the mapper
            evalString = "_mapper = vtk.vtkPolyDataMapper()\n"
            evalString += "_mapper.SetInput(_delaunay.GetOutput())\n"
            evalString += "_mapper.SetLookupTable(_lut)\n"
            # note that zMin and zMax are evaluated in setData()
            evalString += "_mapper.SetScalarRange(_zMin, _zMax)"
            self.renderer.runString(evalString)

        elif self.fname is not None:
            # set up the mapper
            evalString = "_mapper = vtk.vtkDataSetMapper()\n"
            evalString += "_mapper.SetInput(_grid)\n"
            evalString += "_mapper.ScalarVisibilityOn()\n"
            evalString += "_mapper.SetLookupTable(_lut)\n"
            evalString += "_mapper.SetScalarRange(_scalarMin, _scalarMax)"
            self.renderer.runString(evalString)

        # set up the actor
        evalString = "_actor = vtk.vtkActor()\n"
        evalString += "_actor.SetMapper(_mapper)"
        self.renderer.runString(evalString)

        # add the actor to the scene
        self.renderer.runString("_renderer.AddActor(_actor)")

        # set up the text
        # properties (I think this should come from the pyvisi Text() object
        # at some stage, but we'll hard code it here...)
        # I'll also need separate properties for axes, titles, labels etc...
        # but keep them all the same just to get this going
        evalString = "_textProp = vtk.vtkTextProperty()\n"
        evalString += "_textProp.SetFontFamilyToArial()\n"
        evalString += "_textProp.BoldOff()\n"
        evalString += "_textProp.ItalicOff()\n"
        evalString += "_textProp.ShadowOff()\n"
        evalString += "_textProp.SetColor(0,0,0)\n"
        self.renderer.runString(evalString)

        # set the title if set
        if self.title is not None:
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
            evalString += "SetValue(0.5,0.95)\n"# this should be user-settable
            # add the actor to the scene
            evalString += "_renderer.AddActor(_titleActor)"
            self.renderer.runString(evalString)

        # add the axes
        evalString = "_axes = vtk.vtkCubeAxesActor2D()\n"
        evalString += "_axes.SetInput(_grid)\n"
        evalString += "_axes.SetCamera(_renderer.GetActiveCamera())\n"
        evalString += "_axes.SetLabelFormat(\"%6.4g\")\n"
        evalString += "_axes.SetFlyModeToOuterEdges()\n"
        evalString += "_axes.SetFontFactor(0.8)\n"
        evalString += "_axes.SetAxisTitleTextProperty(_textProp)\n"
        evalString += "_axes.SetAxisLabelTextProperty(_textProp)\n"
        evalString += "_axes.GetProperty().SetColor(0,0,0)\n"
        # this next line sets the Z axis visibility off!!  Is a bug in vtk
        # 4.2, dunno how am going to handle this if it is fixed in later
        # versions of vtk
        evalString += "_axes.YAxisVisibilityOff()\n"
        evalString += "_axes.SetNumberOfLabels(5)\n"

        # if we have an xlabel set it
        if self.xlabel is not None:
            evalString += "_axes.SetXLabel(\"%s\")\n" % self.xlabel
        else:
            evalString += "_axes.SetXLabel(\"\")\n"

        # if we have a ylabel set it
        if self.ylabel is not None:
            evalString += "_axes.SetYLabel(\"%s\")\n" % self.ylabel
        else:
            evalString += "_axes.SetXLabel(\"\")\n"

        # add the axes to the scene
        evalString += "_renderer.AddProp(_axes)"
        self.renderer.runString(evalString)

        # add a scalar bar (need to make this an option somewhere!!)
        # I also need to add the ability for the user to set the values of
        # the various parameters set below, and some kind of logical
        # defaults similar to or the same as what I have below.
        evalString = "_scalarBar = vtk.vtkScalarBarActor()\n"
        evalString += "_scalarBar.SetLookupTable(_lut)\n"
        evalString += "_scalarBar.SetWidth(0.1)\n"
        evalString += "_scalarBar.SetHeight(0.7)\n"
        evalString += "_scalarBar.SetPosition(0.9, 0.2)\n"

        # set up the label text properties 
        evalString += "_scalarBarTextProp = _scalarBar.GetLabelTextProperty()\n"
        evalString += "_scalarBarTextProp.ShallowCopy(_textProp)\n"
        evalString += "_scalarBarTextProp.SetFontSize(10)\n"
    
        # add the scalar bar to the scene
        evalString += "_renderer.AddActor(_scalarBar)\n"
        self.renderer.runString(evalString)

        return


# vim: expandtab shiftwidth=4:
