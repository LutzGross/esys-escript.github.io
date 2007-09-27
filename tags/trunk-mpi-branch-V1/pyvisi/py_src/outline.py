"""
@author: John NGUI
"""

import vtk

class Outline:
	"""
	Class that defines an outline.
	"""

	def __init__(self, object):
		"""
		Initialise the outline.

		@type object: vtkUnstructuredGrid, etc
		@param object: Data source to the outline
		"""
		
		self.__object = object
		self.__vtk_outline = vtk.vtkOutlineFilter()
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the outline.
		"""

		self.__vtk_outline.SetInput(self.__object)
		
	
	def _getOutput(self):
		"""
		Return the output of the outline.

		@rtype: vtkPolyData
		@return: Polyognal data
		"""

		return self.__vtk_outline.GetOutput()
	
