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

		# Use arrow instead of cone or sphere.
		self.__vtk_arrow2D.SetGlyphTypeToArrow()
		# Fill the inside of the arrow.
		self.__vtk_arrow2D.SetFilled(0)
	
	def _getOutput(self):
		"""
		Return the output of the 2D arrow.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""
	
		return self.__vtk_arrow2D.GetOutput()


###############################################################################


class Arrow3D:
	"""
	Class that defines 3D arrow.
	"""

	def __init__(self):	
		"""
		Initialise the 3D arrow.
		"""

		self.__vtk_arrow3D = vtk.vtkArrowSource()
		
	def _getOutput(self):
		"""
		Return the output of the 3D arrow.

		@rtype: vtkPolyData
		@return Polygonal data
		"""

		return self.__vtk_arrow3D.GetOutput()

