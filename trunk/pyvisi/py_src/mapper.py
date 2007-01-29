"""
@author: John NGUI
"""

import vtk

class DataSetMapper:
	"""
	Class that defines a data set mapper.
	"""

	# lookup_table = None is used by the outline object.
	def __init__(self, object, lookup_table = None):
		"""
		Initialise the data set mapper.

		@type object: vtkDataSet (i.e. vtkUnstructuredGrid, vtkPolyData, etc) 
		@param object: Data source map
		@type lookup_table: vtkLookupTable
		@param lookup_table: Maps scalar values to a color
		"""

		self.__object = object
		self.__vtk_data_set_mapper = vtk.vtkDataSetMapper()
		self.__setInput()
		if(lookup_table != None):
			self.__setLookupTable(lookup_table)

	def __setInput(self):
		"""
		Set the the input for the data set mapper.
		"""

		self.__vtk_data_set_mapper.SetInput(self.__object)

	def __setLookupTable(self, lookup_table):
		"""
		Set the lookup table for the data set mapper.
	
		@type lookup_table: vtkLookupTable
		@param lookup_table: Map scalar values to a color 
		"""

		self.__vtk_data_set_mapper.SetLookupTable(lookup_table)

	def _setScalarRange(self, range):
		"""
		Set the minimum and maximum scalar range to be mapped into the lookup
		table for the data set mapper.
		
		@type range: tuple
		@param range: Maximum and mimimum scalar range
		"""

		self.__vtk_data_set_mapper.SetScalarRange(range) 

	def _getDataSetMapper(self):
		"""
		Return the data set mapper.

		@rtype: vtkDataSetMapper 	
		@return: Data set mapper
		"""

		return self.__vtk_data_set_mapper

	
