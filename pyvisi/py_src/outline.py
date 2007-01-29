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

		@type object: vtkDataSet (i.e. vtkUnstructuredGrid, vtkPolyData, etc)
		@param object: Data source to outline
		"""
		
		self.__object = object
		self.__vtk_outline = vtk.vtkOutlineFilter()
		self.__setupOutline()

	def __setupOutline(self):
		"""
		Setup the outline.
		"""

		self.__vtk_outline.SetInput(self.__object)
		self.__output = self.__vtk_outline.GetOutput()
		
	
	def _getOutput(self):
		"""
		Return the output of the outline.

		@rtype: vtkPolyData
		@return: Polyognal data
		"""

		return self.__output
	
