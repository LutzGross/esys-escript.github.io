
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"

import vtk
from esys.escript import getMPISizeWorld

class DataSetMapper:
	"""
	Class that defines a data set mapper.
	"""

	def __init__(self):
		"""
		Initialise the data set mapper.
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.DataSetMapper is not running on more than one processor."

		self.__vtk_data_set_mapper = vtk.vtkDataSetMapper()
		# Keeps track whether the scalar range has been specified
		# by the user.
		self.__scalar_range_set = False

	# 'lookup_table = None' is used only by the Outline.
	def _setupDataSetMapper(self, object, lookup_table = None): 
		"""
		Setup the data set mapper.	

		@type object: vtkDataSet (i.e. vtkUnstructuredGrid, vtkPolyData, etc) 
		@param object: Data source map
		@type lookup_table: vtkLookupTable
		@param lookup_table: Maps scalar values to colors
		"""

		self.__object = object
		self.__setInput()

		if(lookup_table != None): # False for the outline.
			self.__setLookupTable(lookup_table)

	def __setInput(self):
		"""
		Set the input for the data set mapper.
		"""

		self.__vtk_data_set_mapper.SetInput(self.__object)

	def __setLookupTable(self, lookup_table):
		"""
		Set the lookup table for the data set mapper.
	
		@type lookup_table: vtkLookupTable
		@param lookup_table: Map scalar values to colors
		"""

		self.__vtk_data_set_mapper.SetLookupTable(lookup_table)

	def setScalarRange(self, lower_range, upper_range):
		"""
		Set the minimum and maximium scalar range for the data set mapper. This
		method is called when the range has been specified by the user. 
		Therefore, the scalar range read from the source will be ignored. 
		
		@type lower_range: Lower range of scalar value
		@param lower_range: Number
		@type upper_range: Upper range of scalar value
		@param upper_range: Number
		"""

		self.__scalar_range_set = True
		self.__vtk_data_set_mapper.SetScalarRange(lower_range, upper_range) 

	def _setScalarRange(self, range):
		"""
		Set the minimum and maximum scalar range for the data set mapper. This
		method is called when the range has NOT been specified by the user. 
		Therefore, the scalar range read from the source will be used instead. 
		
		@type range: Two column tuple containing numbers
		@param range: Minimum and maximum data set mapper scalar range
		"""
		self.__scalar_range_set = True
		self.__vtk_data_set_mapper.SetScalarRange(range) 
	
	def _setScalarVisibilityOn(self):
		"""
		Scalar data is used to color the rendered object.
		"""

		self.__vtk_data_set_mapper.ScalarVisibilityOn()
	
	def _getDataSetMapper(self):
		"""
		Return the data set mapper.

		@rtype: vtkDataSetMapper 	
		@return: Data set mapper
		"""

		return self.__vtk_data_set_mapper
	
	def _getDataSetMapperLookupTable(self):
		"""
		Return the data set mapper's lookup table.

		@rtype: vtkScalarsToColors 	
		@return: Converts scalar data to colors
		"""

		return self.__vtk_data_set_mapper.GetLookupTable()
	
	def _isScalarRangeSet(self):
		"""
		Return whether the data set mapper's scalar range has been specified \
		by the user.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__scalar_range_set 
	
	def _getDataSetMapperRange(self):
		"""
		Return the mapper's scalar range.

		@rtype: Two column tuple containing numbers
		@return: Minimum and maximum range of the data set mapper's scalar range
		"""

		return self.__vtk_data_set_mapper.GetScalarRange()


###############################################################################


class ImageMapper:
	"""
	Class that defines an image mapper.
	"""

	def __init__(self):
		"""
		Initialise the image mapper.
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.ImageMapper is not running on more than one processor."
		self.__vtk_image_mapper = vtk.vtkImageMapper()

	def _setupImageMapper(self, object):
		"""
		Setup the image mapper.
		
		@type object: vtkImageData
		@param object: Image data
		"""

		self.__object = object
		self.__setInput()
		# Both color window and color level needs to be set, otherwise only 
		# a black image will be produced. Both values were obtained from 
		# an example found on the VTK site.
		self.__setColorWindow(255)
		self.__setColorLevel(127.5)

	def __setInput(self):
		"""
		Set the input for the image mapper.
		"""

		self.__vtk_image_mapper.SetInput(self.__object)

	def __setColorWindow(self, color):
		"""
		Set the color of the window.
		"""

		self.__vtk_image_mapper.SetColorWindow(color)

	def __setColorLevel(self, color):
		"""
		Set the color level of the window.
		"""

		self.__vtk_image_mapper.SetColorLevel(color)

	def _getImageMapper(self):
		"""
		Return the image mapper.

		@rtype: vtkImageMapper
		@return: Image mapper
		"""

		return self.__vtk_image_mapper
	
