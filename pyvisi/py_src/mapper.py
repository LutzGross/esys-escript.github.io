"""
@author: John NGUI
"""

import vtk

class DataSetMapper:
	"""
	Class that defines a data set mapper.
	"""

	# 'lookup_table = None' is used only by the Outline.
	def __init__(self):
		"""
		Initialise the data set mapper.
		"""

		self.__vtk_data_set_mapper = vtk.vtkDataSetMapper()

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

	def _setScalarRange(self, range):
		"""
		Set the minimum and maximum scalar range for the data set mapper.
		
		@type range: Two column tuple containing numbers
		@param range: Minimum and maximum data set mapper scalar range
		"""

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


###############################################################################


class ImageMapper:
	"""
	Class that defines a image mapper.
	"""

	def __init__(self):
		"""
		Initialise the image mapper.
		"""

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
