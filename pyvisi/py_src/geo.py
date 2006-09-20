"""
@author: John Ngui
@author: Lutz Gross
"""

class Position:
	"""
	Class that defines the positioning of components in the visualization.
	"""

	def __init__(self, x_coor, y_coor, z_coor):
		"""
		Initialize the x,y and z coordinates.	

		@type x_coor: Number
		@param x_coor: X coordinate in global position
		@type y_coor: Number
		@param y_coor: Y coordinate in global position
		@type z_coor: Number 
		@param z_coor: Z coordinate in global position
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

