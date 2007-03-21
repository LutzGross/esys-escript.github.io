"""
@author: John NGUI
"""

import vtk

class CellDataToPointData:
	"""
	Class that defines a filter to convert cell data to point data by 
	averaging the data values of all cells using a particular point.
	"""

	def __init__(self, object):
		"""
		Initialise the cell to point data filter.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the cell to point data filter
		"""

		self.__object = object
		self.__vtk_cell_to_point = vtk.vtkCellDataToPointData()

		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the cell to point data filter.
		"""

		self.__vtk_cell_to_point.SetInput(self.__object)

	def _getOutput(self):
		"""
		Return the output of the cell to point data filter.

		@rtype: vtkDataSet
		@return: Data set
		"""

		return self.__vtk_cell_to_point.GetOutput()

	
