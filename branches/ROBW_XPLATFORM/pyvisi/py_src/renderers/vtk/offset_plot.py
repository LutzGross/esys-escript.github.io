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

# $Id: offset_plot.py,v 1.1 2005/11/30 03:07:18 paultcochrane Exp $

"""
Class and functions associated with a pyvisi OffsetPlot objects
"""

# generic imports
from common import debugMsg
import numarray
import os
import copy

# module specific imports
from plot import Plot

__revision__ = '$Revision: 1.1 $'

class OffsetPlot(Plot):
    """
    Offset plot
    """
    def __init__(self, scene):
        """
        Initialisation of the OffsetPlot class
        
        @param scene: The Scene to render the plot in
        @type scene: Scene object
        """
        debugMsg("Called OffsetPlot.__init__()")
        Plot.__init__(self, scene)

        self.renderer = scene.renderer
        self.renderer.addToInitStack("# OffsetPlot.__init__()")
        self.renderer.addToInitStack("_plot = vtk.vtkXYPlotActor()")

        self.title = None
        self.xlabel = None
        self.ylabel = None

        # the extra separation between curves (user set)
        self.sep = None

	# the default values for shared info
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

	@param fname: Filename of the input vtk file
	@type fname: string

	@param format: Format of the input vtk file ('vtk' or 'vtk-xml')
	@type format: string

	@param scalars: the name of the scalar data in the vtk file to use
	@type scalars: string
        """
        debugMsg("Called setData() in OffsetPlot()")

        self.renderer.runString("# OffsetPlot.setData()")

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
	    # looks right for an OffsetPlot object

	    raise ImplementationError, "Can't process escript Data yet"

	elif self.otherData:

	    # do some sanity checking on the data
	    if len(dataList) > 3 or len(dataList) < 1:
		raise ValueError, \
			"Must have either one, two or three input arrays"

	    # the data is y values located at different x positions, changing
	    # over time, so the normal x-direction is t, the normal y direction
	    # is both x and y; y basically being offset by the x values
	    # therefore will refer to tData, xData and yData

	    # compare the shapes of the input vectors.
	    # assume that the first one is the t data, and that the first
	    # dimension of the second one is the same length as the t data
	    # length
	    if len(dataList) == 1:
		yData = dataList[0]
	    elif len(dataList) == 2:
		tData = dataList[0]
		yData = dataList[1]
		if tData.shape[0] != yData.shape[0]:
		    raise ValueError, \
			    "Input arrays don't have the correct shape"
	    elif len(dataList) == 3:
		tData = dataList[0]
		xData = dataList[1]
		yData = dataList[2]
		if tData.shape[0] != yData.shape[0]:
		    raise ValueError, \
			"First dim of third arg doesn't agree with first arg"
		if len(yData.shape) == 1:
		    if xData.shape[0] != 1:
			raise ValueError, \
			   "Second arg must be scalar when third arg is vector"
		elif len(yData.shape) == 2:
		    if xData.shape[0] != yData.shape[1]:
			raise ValueError, \
			   "Second dim of 3rd arg doesn't agree with 2nd arg"

	    # if only have one array input, then autogenerate tData
	    if len(dataList) == 1:
		tData = range(1, len(dataList[0])+1)
		if len(tData) != len(dataList[0]):
		    errorString = "Autogenerated xData array length not "
		    errorString += "equal to input array length"
		    raise ValueError, errorString
		## pass around the t data
		self.renderer.renderDict['_t'] = copy.deepcopy(tData)
	    # if have two arrays to plot, the first one is the t data
	    elif len(dataList) == 2:
		tData = dataList[0]
		## pass around the t data
		self.renderer.renderDict['_t'] = copy.deepcopy(tData)
		# don't need the first element of the dataList, so get rid of it
		dataList = dataList[1:]
	    elif len(dataList) == 3:
		## pass around the t data
		self.renderer.renderDict['_t'] = copy.deepcopy(tData)
		## pass around the x data
		self.renderer.renderDict['_x'] = copy.deepcopy(xData)
	    else:
		# shouldn't get to here, but raise an error anyway
		raise ValueError, "Incorrect number of arguments"

	    # set up the vtkDataArray object for the t data
	    self.renderer.runString(
		    "_tData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)")
	    self.renderer.runString(
		    "_tData.SetNumberOfTuples(len(_t))")

	    ## now to handle the y data
	    if len(yData.shape) == 1:
		dataLen = 1
	    elif len(yData.shape) == 2:
		dataLen = yData.shape[1]
	    else:
		raise ValueError, \
			"The last setData argument has the incorrect shape"

	    # share around the y data
	    for i in range(dataLen):
		yDataVar = "_y%d" % i
		if len(yData.shape) == 1:
		    self.renderer.renderDict[yDataVar] = copy.deepcopy(yData)
		else:
		    self.renderer.renderDict[yDataVar] = \
			    copy.deepcopy(yData[:, i])
		# check that the data here is a 1-D array
		if len(self.renderer.renderDict[yDataVar].shape) != 1:
		    raise ValueError, "Can only handle 1D arrays at present"

	elif fname is not None and not self.escriptData and not self.otherData:
	    # now handle the case when we have a file as input
	    raise ImplementationError, "Sorry, can't handle file input yet"

        # concatenate the data
        evalString = "_yAll = concatenate(["
        for i in range(dataLen-1):
            evalString += "_y%d," % i
        evalString += "_y%d])" % int(dataLen-1)
        self.renderer.runString(evalString)

        # grab the min and max values
        self.renderer.runString("_yMax = max(_yAll)")
        self.renderer.runString("_yMin = min(_yAll)")

        # keep the data apart a bit
        if self.sep is None:
            self.renderer.runString("_const = 0.1*(_yMax - _yMin)")
        else:
            evalString = "_const = %f" % self.sep
            self.renderer.runString(evalString)

        # behave differently with the shift if we have xData as to not
        if len(dataList) == 3:
            # this is for when we have xData
            self.renderer.runString("_yMaxAbs = max(abs(_yAll))")
            # calculate the minimum delta x
            x1 = xData[:-1]
            x2 = xData[1:]
            minDeltax = min(x2 - x1)
            evalString = "_scale = %f/(2.0*_yMaxAbs)" % minDeltax
            self.renderer.runString(evalString)

            for i in range(dataLen):
                evalString = "_y%d = _scale*_y%d + _x[%d]" % (i, i, i)
                self.renderer.runString(evalString)
        else:
            # shift the data up
            self.renderer.runString("_shift = _yMax - _yMin + _const")

            for i in range(dataLen):
                evalString = "_y%d = _y%d + %f*_shift" % (i, i, i)
                self.renderer.runString(evalString)

        # set up the vtkDataArray objects
        for i in range(dataLen):
            evalString = \
            "_y%dData = vtk.vtkDataArray.CreateDataArray(vtk.VTK_FLOAT)\n" % i
            evalString += "_y%dData.SetNumberOfTuples(len(_y%d))" % (i, i)
            self.renderer.runString(evalString)

        ## t data
        # put the data into the data arrays
        self.renderer.runString("for _i in range(len(_t)):")
        # need to be careful here to remember to indent the code properly
        evalString = "    _tData.SetTuple1(_i,_t[_i])"
        self.renderer.runString(evalString)

        ## y data
        # put the data into the data arrays
        self.renderer.runString("for _i in range(len(_t)):")
        # need to be careful here to remember to indent the code properly
        for i in range(dataLen):
            evalString = "    _y%dData.SetTuple1(_i,_y%d[_i])" % (i, i)
            self.renderer.runString(evalString)

        for i in range(dataLen):
            # create the field data object
            evalString = "_fieldData%d = vtk.vtkFieldData()\n" % i
            evalString += "_fieldData%d.AllocateArrays(2)\n" % i
            evalString += "_fieldData%d.AddArray(_tData)\n" % i
            evalString += "_fieldData%d.AddArray(_y%dData)" % (i, i)
            self.renderer.runString(evalString)

        for i in range(dataLen):
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
        for i in range(dataLen):
            evalString = "_plot.SetDataObjectYComponent(%d,1)" % i
            self.renderer.runString(evalString)

        # note: am ignoring zlabels as vtk xyPlot doesn't support that
        # dimension for line plots (I'll have to do something a lot more
        # funky if I want that kind of functionality)

        return

    def render(self):
        """
        Does OffsetPlot object specific (pre)rendering stuff
        """
        debugMsg("Called OffsetPlot.render()")

        self.renderer.runString("# OffsetPlot.render()")
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
