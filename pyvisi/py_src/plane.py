"""
@author: John NGUI
"""

import vtk
from position import GlobalPosition

class Plane:
	"""
	Class that defines a plane that cuts/clips rendered objects. 
	"""

	def __init__(self, transform):
		"""
		Initialise the plane.	

		@type transform: L{Transform <transform.Transform>} object
		@param transform: Specifies the orientation of the plane
		"""

		self.__transform = transform
		self.__vtk_plane = vtk.vtkPlane()
		self.__setupPlane()
	
	def __setupPlane(self):	
		"""
		Setup the plane.
		"""

		# Default origin of the of the plane is (0,0,0).
		self.__setOrigin(GlobalPosition(0,0,0))
		# Default normal of the plane is parrallel to the z-axis.
		self.__setNormal(GlobalPosition(0, 0, 1))
		self.__setTransform()

	def __setOrigin(self, position):
		"""
		Set the origin of the plane.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Origin of the plane
		"""
		self.__vtk_plane.SetOrigin(position._getGlobalPosition())

	def __setNormal(self, position):
		"""
		Set the normal of the plane.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Normal of the plane
		"""

		self.__vtk_plane.SetNormal(position._getGlobalPosition())

	def __setTransform(self):
		"""
		Set the transformation of the plane.
		"""

		self.__vtk_plane.SetTransform(self.__transform)

	def _getPlane(self):
		"""
		Return the plane.

		@rtype: vtkPlane
		@return: Plane that cuts through rendered objects
		"""

		return self.__vtk_plane


###############################################################################


class PlaneSource:
	"""
	Class that defines a plane source.
	"""

	def __init__(self):
		"""
		Initialise the plane source. A plane source is defined by an origin
		and two other points, which form the axes. By default the origin is 
		set to (0, 0, 0) and the two points are set to (1, 0, 0) and (0, 1, 0).
		"""

		self.__vtk_plane_source = vtk.vtkPlaneSource()
		# Default origin of the plane source is (0, 0, 0).
		self.__setOrigin(GlobalPosition(0, 0, 0))
		# Default position of the two point from the origin is (1, 0, 0) and 
		# (0, 1, 0).
		self.setPoint1(GlobalPosition(1, 0, 0))
		self.setPoint2(GlobalPosition(0, 1, 0))

	def __setOrigin(self, position):
		"""
		Set the origin of the plane source.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Plane source origin
		"""

		self.__vtk_plane_source.SetOrigin(position._getGlobalPosition())

	def setPoint1(self, position):
		"""
		Set the first point from the origin of the plane source.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Position of the first axis of the plane source
		"""

		self.__vtk_plane_source.SetPoint1(position._getGlobalPosition())

	def setPoint2(self, position):
		"""
		Set the second point from the origin of the plane source.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Position of the second axis of the plane source
		"""

		self.__vtk_plane_source.SetPoint2(position._getGlobalPosition())


	def _getOutput(self):
		"""
		Return the output of the plane source.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_plane_source.GetOutput()

