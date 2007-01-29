"""
@author: John NGUI
"""

import vtk

class Normals:
	"""
	Class that defines normals. Normals is used to average the normals of 
	points in order to generate better sufaces (in the case of tensors, normals
	avoids the tensors from appearing black in color).
	"""

	def __init__(self, object):
		"""
		Initialise the normals.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the normals
		"""

		self.__object = object
		self.__vtk_poly_data_normals = vtk.vtkPolyDataNormals()
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the normals.
		"""

		self.__vtk_poly_data_normals.SetInput(self.__object)

	def _getOutput(self):
		"""
		Return the output of the normals.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_poly_data_normals.GetOutput()
		
