"""
@author: John NGUI
"""

import vtk

class Normals:
	"""
	Class that defines normals. Normals are used to average the normals of 
	points in order to generate better sufaces (in the case of tensors, normals
	avoids the tensors from appearing black in color).
	"""

	def __init__(self):
		"""
		Initialise the normals.
		"""

		self.__vtk_poly_data_normals = vtk.vtkPolyDataNormals()

	def _setupNormals(self, object):
		"""
		Setup the normals.

		@type object: vtkPolyData, etc
		@param object: Input for the normals
		"""

		self.__object = object
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the normals.
		"""

		self.__vtk_poly_data_normals.SetInput(self.__object)

	def _getNormalsOutput(self):
		"""
		Return the output of the normals.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_poly_data_normals.GetOutput()
