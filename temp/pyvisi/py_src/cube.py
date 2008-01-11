"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


import vtk

class CubeSource:
	"""
	Class that defines a cube source. The center of the cube souce defines
	the point from which the cube is to be generated and the X, Y 
	and Z lengths define the length of the cube from the center point. If
	X length is 3, then the X length to the left and right of the center 
	point is 1.5 respectively.
	"""
	
	def __init__(self):
		"""
		Initialise the cube source. 
		"""

		self.__vtk_cube_source = vtk.vtkCubeSource()

	def setCenter(self, center):
		"""
		Set the cube source center.

		@type center: L{GlobalPosition <position.GlobalPosition>} object
		@param center: Cube source center
		"""

		self.__vtk_cube_source.SetCenter(center._getGlobalPosition())

	def setXLength(self, length):
		"""
		Set the cube source length along the x-axis.

		@type length: Number
		@param length: Cube source length along the x-axis
		"""
		self.__vtk_cube_source.SetXLength(length)

	def setYLength(self, length):
		"""
		Set the cube source length along the y-axis.

		@type length: Number
		@param length: Cube source length along the y-axis
		"""

		self.__vtk_cube_source.SetYLength(length)

	def setZLength(self, length):
		"""
		Set the cube source length along the z-axis.

		@type length: Number
		@param length: Cube source length along the z-axis
		"""

		self.__vtk_cube_source.SetZLength(length)

	def _getCubeSourceOutput(self):
		"""
		Return the output of the cube source

		@rtype: vtkPolyData
		@return: Polygonal data
		"""
	
		return self.__vtk_cube_source.GetOutput()



