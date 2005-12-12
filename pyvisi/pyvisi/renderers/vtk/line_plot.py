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

# $Id: line_plot.py,v 1.1 2005/11/30 03:07:18 paultcochrane Exp $

"""
Class and functions associated with a pyvisi LinePlot objects
"""

# generic imports
from pyvisi.renderers.vtk.common import debugMsg
import Numeric
import os
import copy

# module specific imports
from pyvisi.renderers.vtk.plot import Plot

__revision__ = '$Revision: 1.1 $'

class LinePlot(Plot):
    """
    Line plot
    """
    def __init__(self, scene):
        """
        Initialisation of the LinePlot class
        
        @param scene: The Scene to render the plot in
        @type scene: Scene object
        """
        debugMsg("Called LinePlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer
        self.renderer.addToInitStack("# LinePlot.__init__()")
        self.renderer.addToInitStack("_plot = vtk.vtkXYPlotActor()")

        # offset the data in the lineplot?
        self.offset = False

	# default values for stuff to be passed around
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

	@param offset: whether or not to offset the lines from one another
	@type offset: boolean

	@param fname: Filename of the input vtk file
	@type fname: string

	@param format: format of the input vtk file ('vtk' or 'vtk-xml')
	@type format: string

	@param scalars: the name of the scalar data in the vtk file to use
	@type scalars: string
        """
        debugMsg("Called setData() in LinePlot()")

        self.renderer.runString("# LinePlot.setData()")

	# process the options, if any
	## offset
        if options.has_key('offset'):
            self.offset = options['offset']
        else:
            self.offset = False
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
		escriptY = dataList[0]
		escriptX = escriptY.getDomain().getX()
	    elif len(dataList) == 2:
		# ok, have two data variables, first is x data
		escriptX = dataList[0]
		escriptY = dataList[1]
	    else:
		errorString = \
			"Expecting 1 or 2 elements in data list.  I got %d" \
			% len(dataList)
		raise ValueError, errorString

	    ####!!!! now process the data properly so that I can plot it
	    print "escriptX.shape() = %s" % escriptX.shape()
	    print "escriptY.shape() = %s" % escriptY.shape()

	    # need to check the length and shape of the data to make sure it
	    # looks right for a LinePlot object

	    raise ImplementationError, "Can't process escript Data yet"

	elif self.otherData:
	    # do some sanity checking on the data
	    for i in range(len(dataList)):
		if len(dataList[0]) != len(dataList[i]):
		    raise ValueError, \
			    "Input vectors must all be the same length"

	    # if have more than one array to plot, the first one is the x data
	    if len(dataList) > 1:
		xData = dataList[0]
		## pass the x data around
		self.renderer.renderDict['_x'] = copy.deepcopy(xData)
		# don't need the first element of the dataList, so get rid of it
		dataList = dataList[1:]
		# if only have one array input, then autogenerate xData
	    elif len(dataList) == 1:
		xData = range(1, len(dataList[0])+1)
		if len(xData) != len(dataList[0]):
		    errorString = "Autogenerated xData array length not "
		    errorString += "equal to input array length"
		    raise ValueError, errorString
		## pass the x data around
		self.renderer.renderDict['_x'] = copy.deepcopy(xData)

	    # set up the vtkDataArray object for the x data
	    self.renderer.runString(
		    "_xData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)")
	    self.renderer.runString(
		    "_xData.SetNumberOfTuples(len(_x))")

	    ## now to handle the y data

	    # now to add my dodgy hack until I have a decent way of sharing data
	    # objects around properly
	    for i in range(len(dataList)):
		# check that the data here is a 1-D array
		if len(dataList[i].shape) != 1:
		    raise ValueError, "Can only handle 1D arrays at present"

		yDataVar = "_y%d" % i
		self.renderer.renderDict[yDataVar] = copy.deepcopy(dataList[i])

	elif fname is not None and not self.escriptData and not self.otherData:
	    # now handle the case when we have a file as input
	    raise ImplementationError, "Sorry, can't handle file input yet"

        # if offset is true then shift the data
        if self.offset:
            # concatenate the data
            evalString = "_yAll = concatenate(["
            for i in range(len(dataList)-1):
                evalString += "_y%d," % i
            evalString += "_y%d])" % int(len(dataList)-1)
            self.renderer.runString(evalString)

            # grab the min and max values
            self.renderer.runString("_yMax = max(_yAll)")
            self.renderer.runString("_yMin = min(_yAll)")

            # keep the data apart a bit
            self.renderer.runString("_const = 0.1*(_yMax - _yMin)")

            # now shift the data
            self.renderer.runString("_shift = _yMax - _yMin + _const")
            for i in range(len(dataList)):
                evalString = "_y%d = _y%d + %d*_shift" % (i, i, i)
                self.renderer.runString(evalString)

        # set up the vtkDataArray objects
        for i in range(len(dataList)):
            evalString = \
            "_y%dData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)\n" % i
            evalString += "_y%dData.SetNumberOfTuples(len(_y%d))" % (i, i)
            self.renderer.runString(evalString)

        ## x data
        # put the data into the data arrays
        self.renderer.runString("for _i in range(len(_x)):")
        # need to be careful here to remember to indent the code properly
        evalString = "    _xData.SetTuple1(_i,_x[_i])"
        self.renderer.runString(evalString)

        ## y data
        # put the data into the data arrays
        self.renderer.runString("for _i in range(len(_x)):")
        # need to be careful here to remember to indent the code properly
        for i in range(len(dataList)):
            evalString = "    _y%dData.SetTuple1(_i,_y%d[_i])" % (i, i)
            self.renderer.runString(evalString)

        for i in range(len(dataList)):
            # create the field data object
            evalString = "_fieldData%d = vtk.vtkFieldData()" % i
            self.renderer.runString(evalString)
            evalString = "_fieldData%d.AllocateArrays(2)" % i
            self.renderer.runString(evalString)
            evalString = "_fieldData%d.AddArray(_xData)" % i
            self.renderer.runString(evalString)
            evalString = "_fieldData%d.AddArray(_y%dData)" % (i, i)
            self.renderer.runString(evalString)

        for i in range(len(dataList)):
            # now put the field data into a data object
            evalString = "_dataObject%d = vtk.vtkDataObject()\n" % i
            evalString += "_dataObject%d.SetFieldData(_fieldData%d)\n" % (i, i)

            # the actor should be set up, so add the data object to the actor
            evalString += "_plot.AddDataObjectInput(_dataObject%d)" % i
            self.renderer.runString(evalString)

        # tell the actor to use the x values for the x values (rather than
        # the index)
        self.renderer.runString("_plot.SetXValuesToValue()")

        # set which parts of the data object are to be used for which axis
        self.renderer.runString("_plot.SetDataObjectXComponent(0,0)")
        for i in range(len(dataList)):
            evalString = "_plot.SetDataObjectYComponent(%d,1)" % i
            self.renderer.runString(evalString)

        # note: am ignoring zlabels as vtk xyPlot doesn't support that
        # dimension for line plots (I'll have to do something a lot more
        # funky if I want that kind of functionality)

        # should this be here or elsewhere?
        evalString = "_plot.GetXAxisActor2D().GetProperty().SetColor(0, 0, 0)\n"
        evalString += \
		"_plot.GetYAxisActor2D().GetProperty().SetColor(0, 0, 0)\n"
        evalString += "_renderer.SetBackground(1.0, 1.0, 1.0)"
        self.renderer.runString(evalString)

        # set up the lookup table for the appropriate range of colours
        evalString = "_lut = vtk.vtkLookupTable()\n"
        evalString += "_lut.Build()\n"
        evalString += "_colours = []\n"
        # need to handle the case when only have one element in dataList
        if len(dataList) == 1:
            evalString += "_colours.append(_lut.GetColor(0))\n"
        else:
            for i in range(len(dataList)):
                evalString += "_colours.append(_lut.GetColor(%f))\n" \
                        % (float(i)/float(len(dataList)-1),)
        self.renderer.runString(evalString)
    
        # change the colour of the separate lines
        for i in range(len(dataList)):
            evalString = "_plot.SetPlotColor(%d, _colours[%d][0], " % (i, i)
            evalString += "_colours[%d][1], _colours[%d][2])" % (i, i)
            self.renderer.runString(evalString)

        # make sure the plot is a decent size
        # the size of the actor should be 80% of the render window
        evalString = "_plot.SetPosition(0.1, 0.1)\n" # (0.1 = (1.0 - 0.8)/2)
        evalString += "_plot.SetWidth(0.8)\n"
        evalString += "_plot.SetHeight(0.8)"
        self.renderer.runString(evalString)

        return

    def render(self):
        """
        Does LinePlot object specific (pre)rendering stuff
        """
        debugMsg("Called LinePlot.render()")

        self.renderer.runString("# LinePlot.render()")
        self.renderer.runString("_renderer.AddActor2D(_plot)")

        # set the title if set
        if self.title is not None:
            evalString = "_plot.SetTitle(\'%s\')" % self.title
            self.renderer.runString(evalString)

        # if an xlabel is set, add it
        if self.xlabel is not None:
            evalString = "_plot.SetXTitle(\'%s\')" % self.xlabel
            self.renderer.runString(evalString)

        # if an ylabel is set, add it
        if self.ylabel is not None:
            evalString = "_plot.SetYTitle(\'%s\')" % self.ylabel
            self.renderer.runString(evalString)

        return


# vim: expandtab shiftwidth=4:
