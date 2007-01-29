"""
@author: John NGUI
"""

import vtk

class Transform:
	"""
	Class that defines the orientation of rendered objects.
	"""

	def __init__(self):
		"""
		Initialise the transform object.
		"""

		# NOTE: VTK's values are not accurate. Origin is not exaclty (0,0,0)
		# and normal is not exactly (0, 0, 1). There is a slight 
		# variance. As a result, a slight alteration has to be done in order 
		# for the rendered object to be displayed correctly. Otherwise, the 
		# rendered object may just fall outside the bounding box and nothing 
		# is displayed.  

		self.__OFFSET_VARIANCE = 0.0000000001
		self.__vtk_transform = vtk.vtkTransform()

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

		self.__vtk_transform.Translate(-x_offset, -y_offset, -z_offset)
	
	def rotateX(self, angle):
		"""
		Rotate the rendered object along the x-axis.
		@type angle: Number
		@param angle: Angle to rotate the camera
		"""

		self.__vtk_transform.RotateX(-angle)

	def rotateY(self, angle):
		"""
		Rotate the rendered object along the y-axis.
		@type angle: Number
		@param angle: Angle to rotate the camera
		"""

		self.__vtk_transform.RotateY(angle)
				

	def rotateZ(self, angle):
		"""
		Rotate the rendered object along the z-axis.
		@type angle: Number
		@param angle: Angle to rotate the camera
		"""

		self.__vtk_transform.RotateZ(angle)

	def setPlaneToXY(self, offset = 0):
		"""
		Set the plane orthogonal to the z-axis.
		@type offset: Number
		@param offset: Amount to translate
		"""
	
		self.translate(0, 0, offset + self.__OFFSET_VARIANCE)

	def setPlaneToYZ(self, offset = 0):
		"""
		Set the plane orthogonal to the x-axis.
		@type offset: Number
		@param offset: Amount to translate
		"""
		
		# NOTE: rotateY must come first before translate. Otherwise, 
		# the output may be incorrect.
		self.rotateY(90) 
		self.translate(offset, 0, 0)

	def setPlaneToXZ(self, offset = 0):
		"""
		Set the plane orthogonal to the y-axis.
		@type offset: Number
		@param offset: Amount to translate
		"""

		# rotateX must come first before translate. Otherwise, it won't work.
		self.rotateX(90)
		self.translate(0, offset, 0)
	
	def _getTransform(self):
		"""
		Return the transform instance.
		@rtype: vtkTransform
		@return: VTK transform that is used to specify the orientation
			of objects
		"""

		return self.__vtk_transform
