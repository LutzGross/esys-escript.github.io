"""
@author: John NGUI
"""

import vtk

class Sphere:
	"""
	Class that defines a sphere.
	"""

	def __init__(self):
		"""
		Initialise the sphere.
		"""

		self.__vtk_sphere = vtk.vtkSphereSource()

	def setThetaResolution(self, resolution):
		"""
		Set the theta resolution of the sphere.

		@type resolution: Number
		@param resolution: Theta resolution
		"""

		self.__vtk_sphere.SetThetaResolution(resolution)

	def setPhiResolution(self, resolution):
		"""
		Set the phi resolution of the sphere.

		@type resolution: Number
		@param resolution: Phi resolution
		"""

		self.__vtk_sphere.SetPhiResolution(resolution)

	def _getOutput(self):
		"""
		Return the output of the sphere.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_sphere.GetOutput()

