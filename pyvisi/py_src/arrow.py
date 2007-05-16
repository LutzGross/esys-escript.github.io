"""
@author: John NGUI
"""

import vtk

class Arrow2D:
	"""
	Class that defines 2D arrows.
	"""
	
	def __init__(self):
		"""
		Initialise the 2D arrows.
		"""

		self.__vtk_arrow2D = vtk.vtkGlyphSource2D()
		self.__setupArrow2D()

	def __setupArrow2D(self):
		"""
		Setup the 2D arrows.
		"""

		# Use arrows instead of cone or sphere.
		self.__vtk_arrow2D.SetGlyphTypeToArrow()
		# Fill the inside of the arrows.
		self.__vtk_arrow2D.SetFilled(0)
	
	def _getArrow2DOutput(self):
		"""
		Return the output of the 2D arrows.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""
	
		return self.__vtk_arrow2D.GetOutput()


###############################################################################


class Arrow3D:
	"""
	Class that defines 3D arrows.
	"""

	def __init__(self):	
		"""
		Initialise the 3D arrows.
		"""

		self.__vtk_arrow3D = vtk.vtkArrowSource()
		
	def _getArrow3DOutput(self):
		"""
		Return the output of the 3D arrows.

		@rtype: vtkPolyData
		@return Polygonal data
		"""

		return self.__vtk_arrow3D.GetOutput()

