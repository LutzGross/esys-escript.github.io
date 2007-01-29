"""
@author: John NGUI
"""

import vtk

class Arrow2D:
	"""
	Class that defines 2D arrow.
	"""
	
	def __init__(self):
		"""
		Initialise the 2D arrow.
		"""

		self.__vtk_arrow2D = vtk.vtkGlyphSource2D()
		self.__setupArrow2D()

	def __setupArrow2D(self):
		"""
		Setup the 2D arrow.
		"""

		# Set to use arrows instead of cone, sphere, etc.
		self.__vtk_arrow2D.SetGlyphTypeToArrow()
		# Fill the inside of the arrow.
		self.__vtk_arrow2D.SetFilled(0)
		self.__output = self.__vtk_arrow2D.GetOutput()
	
	def _getOutput(self):
		"""
		Return the 2D arrow.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""
	
		return self.__output 

class Arrow3D:
	"""
	Class that defines 3D arrow.
	"""

	def __init__(self):	
		"""
		Initialise the 3D arrow.
		"""

		self.__vtk_arrow3D = vtk.vtkArrowSource()
		self.__output = self.__vtk_arrow3D.GetOutput()
		
	def _getOutput(self):
		"""
		Return the 3D arrow.

		@rtype: vtkPolyData
		@return Polygonal data
		"""

		return self.__output
		

