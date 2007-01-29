"""
@author: John NGUI
"""

import vtk
from position import GlobalPosition

class Plane:
	"""
	Class that defines a plane that cuts through meshes. 
	"""

	def __init__(self, transform):
		"""
		Initialise the plane.	

		@type transform: L{Transform <transform.Transform>} object
		@param transform: Specifies the orientation of the rendered plane
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
		#self.__setOrigin(GlobalPosition(center[0], center[1], center[2]))
		# Default normal of the plane is parrallel to the z-axis.
		self.__setNormal(GlobalPosition(0, 0, 1))
		#self.__setNormal(GlobalPosition(center[0], center[1], center[2]*1000))
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
		Set the transformation to the plane.
		"""

		self.__vtk_plane.SetTransform(self.__transform)

	def _getPlane(self):
		"""
		Return the plane.

		@rtype: vtkPlane
		@return: Plane that cuts through meshes
		"""

		return self.__vtk_plane


class PlaneSource:
	"""
	Class that defines a plane source.
	"""

	def __init__(self):
		"""
		Initialise the plane source.
		"""

		self.__vtk_plane_source = vtk.vtkPlaneSource()
	
	def _getOutput(self):
		"""
		Return the plance source.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_plane_source.GetOutput()


