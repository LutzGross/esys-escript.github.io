"""
Class that defines the positioning of components in the visualization.
"""

class Position:
	"""
	@author: John Ngui
	@author: Lutz Gross
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



#def Position(object):
"""
A position in global coordinates
"""
pass

def Origin(Position):
    """
    The position of the origin
    """
    pass

def Direction(object):
    """
    A dirction in global coordinates
    """
    pass

def XAxis(Direction):
    """
    The direction of the x-axis
    """
    pass

def YAxis(Direction):
    """
    The direction of the y-axis
    """
    pass

def ZAxis(Direction):
    """
    The direction of the z-axis
    """
    pass

def Plane(object):
    """
    A plane in global coordinates
    """
    pass

def XYPlane(Plane):
    """
    The XY plane orthogonal to the z-axis
    """
    pass

def YZPlane(Plane):
    """
    The YZ plane orthogonal to the x-axis
    """
    pass

def ZXPlane(Plane):
    """
    The ZX plane orthogonal to the y-axis
    """
    pass

def Sphere(object):
    """
    A sphere
    """
    pass

