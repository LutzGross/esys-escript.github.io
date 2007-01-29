"""
@author: John NGUI
"""

import vtk

class Cutter:
	"""
	Class that defines a cutter.
	"""

	def __init__(self, object, plane):
		"""
		Initialise the cutter.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the cutter
		@type plane: vtkPlane
		@param plane: Plane to cut the object
		"""

		self.__object = object
		self.__plane = plane
		self.__vtk_cutter = vtk.vtkCutter()

		self.__setInput()
		self.__setCutFunction()
		self.__vtk_cutter.Update()

	def __setInput(self):
		"""
		Set the input for the cutter.
		"""

		self.__vtk_cutter.SetInput(self.__object)

	def __setCutFunction(self):
		"""
		Set the cut functions.
		"""

		self.__vtk_cutter.SetCutFunction(self.__plane)

	def _getOutput(self):
		"""
		Return the output of the cutter.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_cutter.GetOutput()
