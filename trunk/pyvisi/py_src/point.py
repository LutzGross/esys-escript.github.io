"""
@author: John NGUI
"""

import vtk
from position import GlobalPosition

class PointSource:
	"""
	Class that defines the source (location) to generate points. The points are
	generated within the radius of a sphere.
	"""

	def __init__(self, object):
		"""
		Initialise the point source.
		"""

		self.__object = object
		self.__vtk_point_source = vtk.vtkPointSource()

		self.__setupPointSource()

	def __setupPointSource(self):
		"""
		Setup the point source.
		"""

		# Default number of points to generate is 10.
		self.setNumberOfPoints(10)
		# Default center is of the sphere is the center of the data object.
		center = self.__object.GetCenter()
		self.setCenter(GlobalPosition(center[0], center[1], center[2]))
		# Default radius of the sphere is 0.1.
		self.setPointSourceRadius(0.1)
		self.__vtk_point_source.Update()

	def setPointSourceRadius(self, radius):
		"""
		Set the radius of the sphere.

		@type radius: Number
		@param radius: Radius of the sphere
		"""

		self.__vtk_point_source.SetRadius(radius)

	def setCenter(self, position):
		"""
		Set the center of the sphere.
		
		@type position: L{GLobalPosition <position.GlobalPosition>} object
		@param position: Center of the sphere
		"""

		self.__vtk_point_source.SetCenter(position._getGlobalPosition())

	def setNumberOfPoints(self, points):
		"""
		Set the number of points to generate within the sphere.

		@type points: Number
		@param points: Number of points to generate
		"""

		self.__vtk_point_source.SetNumberOfPoints(points)	

	def _getOutput(self):
		"""
		Return the point source.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_point_source.GetOutput()


class StructuredPoints:
	"""
	Class that defines the structured points.
	"""
	
	# NOTE: The algorithm of this class was extracted from Mayavi's 
	# online source code.
	def __init__(self, object):
		"""
		Initialise the structured points.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the structured points
		"""

		self.__object = object
		self.__vtk_structured_points = vtk.vtkStructuredPoints()
		self.__setupStructuredPoints()

	def __setupStructuredPoints(self):
		"""
		Setup the structured points.
		"""

		b = self.__object.GetBounds()
		self.__vtk_structured_points.SetOrigin(b[0], b[2], b[4])
		self.__l = [b[1] - b[0], b[3] - b[2], b[5] - b[4]]
		tot_len = float(self.__l[0] + self.__l[1] + self.__l[2])
		npnt = pow(self.__object.GetNumberOfPoints(), 1. / 3.) + 0.5
		self.__fac = 3.0*npnt / tot_len

		# Default dimension is 1, 1, 1.
		self.setDimension(1, 1, 1)

	def __setExtent(self, x0, x1, y0, y1, z0, z1):
		"""
		Set the extent in the order of x, y, z axes. 

		@type x0: Index of the first point on the x-axis
		@param x1: Index of the last point on the x-axis
		@type y0: Index of the first point on the y-axis
		@param y1: Index of the last point on the y-axis
		@type z0: Index of the first point on the z-axis
		@param z1: Index of the last point on the z-axis
		"""

		self.__vtk_structured_points.SetExtent(x0, x1, y0, y1, z0, z1)

	def __setUpdateExtent(self, x0, x1, y0, y1, z0, z1):
		"""
		Set the update extent in the oder of x, y, z axes.

		@type x0: Number
		@param x1: Index of the last point on the x-axis
		@type y0: Number
		@param y1: Index of the last point on the y-axis
		@type z0: Number
		@param z1: Index of the last point on the z-axis
		"""

		self.__vtk_structured_points.SetUpdateExtent(x0, x1, y0, y1, z0, z1)

	#def __setWholeExtent(self, x0, x1, y0, y1, z0, z1):
	#	self.__vtk_structured_points.SetWholeExtent(x0, x1, y0, y1, z0, z1)

	def setDimension(self, x, y, z):
		"""
		Set the dimension (number of points) on the x, y and z axes. The
		smaller the dimension, the more points are populated.

		@type x: Number
		@param x: Number of points on the x-axis
		@type y: Number
		@param y: Number of points on the y-axis
		@type z: Number
		@param z: Number of points on the z-axis
		"""

		self.__dims = [int((self.__l[0]*self.__fac)/x) + 1, 
				int((self.__l[1] * self.__fac) / y) + 1, 
				int((self.__l[2] * self.__fac) / z) + 1]

		self.__setExtent(0, self.__dims[0] - 1, 0, self.__dims[1] - 1, 0, 
				self.__dims[2] - 1)
		self.__setUpdateExtent(0, self.__dims[0] - 1, 0, self.__dims[1] - 1, 0,
				self.__dims[2] - 1)

		#self.__setWholeExtent(0, self.__dims[0] - 1, 0, self.__dims[1] - 1, 0, 
		#		self.__dims[2] - 1)

		self.__vtk_structured_points.SetDimensions(self.__dims)
		self.__setSpacing()

	def __setSpacing(self):
		"""
		Set the spacing (width, height and length) of the cells that make up
		the dataset.
		"""

		self.__dims = [max(1, x - 1) for x in self.__dims]
		self.__l = [max(1e-3, x) for x in self.__l]
		sp = [self.__l[0] / self.__dims[0], self.__l[1] / self.__dims[1], 
				self.__l[2] / self.__dims[2]]

		self.__vtk_structured_points.SetSpacing(sp)
		# Update the changes made.
		self.__vtk_structured_points.Update()

	def _getStructuredPoints(self):
		"""
		Return the structured points.

		@rtype: vtkStructuredPoints
		@return: Structured points
		"""

		return self.__vtk_structured_points

