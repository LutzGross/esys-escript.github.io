
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


class LocalPosition:
	"""
	Class that defines the local positioning (X and Y)coordiante system (2D).
	"""

	def __init__(self, x_coor, y_coor):
		"""
		Initialise the local position.

		:type x_coor: Number
		:param x_coor: x coordinate
		:type y_coor: Number
		:param y_coor: y coordinate
		"""

		self.__x_coor = x_coor
		self.__y_coor = y_coor
		self.__position = [x_coor, y_coor]

	def _getXCoor(self):
		"""
		Return the X coordinate.

		:rtype: Number
		:return: X coordinate
		"""

		return self.__x_coor

	def _getYCoor(self):
		"""
		Return the Y coordinate.

		:rtype: Number
		:return: Y coordinate
		"""

		return self.__y_coor

	def _getLocalPosition(self):
		"""
		Return the local position.

		:rtype: Two column tuple containing numbers
		:return: Tuple with the x and y coordinates
		"""

		return self.__position


###############################################################################


class GlobalPosition:
	"""
	Class that defines the global positioning (X, Y and Z) coordinate system 
	(3D)
	"""
	
	def __init__(self, x_coor, y_coor, z_coor):
		"""
		Initialise the global position.

		:type x_coor: Number
		:param x_coor: x coordinate
		:type y_coor: Number
		:param y_coor: y coordinate
		:type z_coor: Number
		:param z_coor: z coordinate
		"""

		self.__position = [x_coor, y_coor, z_coor]

	def _getGlobalPosition(self):
		"""
		Return the global position.

		:rtype: Three column tuple containing numbers
		:return: Tuple with the x, y and z coordinates
		"""

		return self.__position
