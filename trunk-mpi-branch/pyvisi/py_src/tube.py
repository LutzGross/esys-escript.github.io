"""
@author: John NGUI
"""

import vtk

class Tube:
	"""
	Class that defines tubes that are wrapped around the streamlines.
	"""

	def __init__(self, object):
		"""
		Initialise the tube.
		"""

		self.__object = object
		self.__vtk_tube = vtk.vtkTubeFilter()

		self.__setupTube()

	def __setupTube(self):
		"""
		Setup the tube.
		"""

		self.__setInput()
		# Default radius of the tube is 0.02.
		self.setTubeRadius(0.02)
		# Default number of sides for the tube is 12.
		self.setTubeNumberOfSides(12)
		#self.setRadiusToVaryByVector()
		self.__vtk_tube.Update()

	def __setInput(self):
		"""
		Set the input for the tube.
		"""

		self.__vtk_tube.SetInput(self.__object)

	def setTubeRadius(self, radius):
		"""
		Set the radius of the tube.
		
		@type radius: Number
		@param radius: Radius of the tube
		"""

		self.__vtk_tube.SetRadius(radius)

	def setTubeNumberOfSides(self, sides):
		"""
		Set the number of sides for the tube. Minimum number of sides is 3.
		The larger the number of sides, the higher quality the steramlines are.
		
		@type sides: Number
		@param sides: Number of sides for the tube
		"""

		self.__vtk_tube.SetNumberOfSides(sides)

	def setTubeRadiusToVaryByVector(self):
		"""
		Set the radius to vary by vector data.
		"""

		self.__vtk_tube.SetVaryRadiusToVaryRadiusByVector()
	
	def setTubeRadiusToVaryByScalar(self):
		"""
		Set the radius to vary by scalar data.
		"""

		self.__vtk_tube.SetVaryRadiusToVaryRadiusByScalar()

	def _getOutput(self):
		"""
		Return the output of the tube.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_tube.GetOutput()
