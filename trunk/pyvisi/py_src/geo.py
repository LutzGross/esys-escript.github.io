"""
@author: John Ngui
@author: Lutz Gross
"""

class Position:
	"""
	Class that defines the x, y and z coordinates of components.
	"""

	def __init__(self, x_coor, y_coor, z_coor):
		"""
		@type x_coor: Number
		@param x_coor: X coordinate 
		@type y_coor: Number
		@param y_coor: Y coordinate
		@type z_coor: Number 
		@param z_coor: Z coordinate
		"""
		
		self.x_coor = x_coor 
		self.y_coor = y_coor 
		self.z_coor = z_coor

	def getXCoor(self):
		"""	
		Return the x coordinate.

		@rtype: Number
		@return: X coordinate
		"""
		return self.x_coor

	def getYCoor(self):
		"""
		Return the y coordinate.

		@rtype: Number
		@return: Y coordiante
		"""

		return self.y_coor

	def getZCoor(self):
		"""
		Return the z coordinate

		@rtype: Number
		@return: Z coordinate
		"""

		return self.z_coor

import vtk

class Transform:
	"""
	Class that defines the orientation of rendered objects.
	"""

	def __init__(self):
		self.vtk_transform = vtk.vtkTransform()
		# Set the transformation to occur after any transformations 
		# represented by the current matrix.	
		#self.vtk_transform.PostMultiply()

	def translate(self, x_offset, y_offset, z_offset):
		"""
		Translate the rendered object along the x, y and z-axes.
		@type x_offset: Number
		@param x_offset: Amount to translate along the x-axis
		@type y_offset: Number
		@param y_offset: Amount to translate along the y-axis
		@type z_offset: Number
		@param z_offset: Amount to translate along the z-axis 
		"""

		self.vtk_transform.Translate(-x_offset, -y_offset, -z_offset)
	
	def rotateX(self, angle):
		"""
		Rotate the rendered object along the x-axis.
		@type angle: Number
		@param angle: Angle to rotate the camera
		"""

		self.vtk_transform.RotateX(-angle)

	def rotateY(self, angle):
		"""
		Rotate the rendered object along the y-axis.
		@type angle: Number
		@param angle: Angle to rotate the camera
		"""

		self.vtk_transform.RotateY(angle)

	def rotateZ(self, angle):
		"""
		Rotate the rendered object along the z-axis.
		@type angle: Number
		@param angle: Angle to rotate the camera
		"""

		self.vtk_transform.RotateZ(angle)

	def xyPlane(self, offset = 0):
		"""
		Set the plane orthogonal to the z-axis.
		@type offset: Number
		@param offset: Amount to translate
		"""

		self.translate(0, 0, offset)

	def yzPlane(self, offset = 0):
		"""
		Set the plane orthogonal to the x-axis.
		@type offset: Number
		@param offset: Amount to translate
		"""
		
		# rotateY must come first before translate. Otherwise, it won't work.
		self.rotateY(89.9) 
		self.translate(offset, 0, 0)

	def xzPlane(self, offset = 0):
		"""
		Set the plane orthogonal to the y-axis.
		@type offset: Number
		@param offset: Amount to translate
		"""

		# rotateX must come first before translate. Otherwise, it won't work.
		self.rotateX(89.9)
		self.translate(0, offset, 0)
	
	def getTransform(self):
		"""
		Return the transform instance.
		@rtype: vtkTransform
		@return: VTK transform that is used to specify the orientation
			of objects
		"""

		return self.vtk_transform
