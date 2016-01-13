
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
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class ContourModule:
	"""
	Class that defines the contour module.
	"""

	def __init__(self):
		"""
		Initliase the contour module.
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.ContourModule is not running on more than one processor."
		self.__vtk_contour = vtk.vtkContourFilter()
		# Keeps track whether the number of contours and its range have 
		# been specified.
		self.__contours = None
		self.__lower_range = None
		self.__upper_range = None

	def _setupContourModule(self, object):
		"""
		Setup the contour module.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the contour
		"""

		self.__object = object
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the contour.
		"""

		self.__vtk_contour.SetInput(self.__object)

	# This method is used to delay the execution of generating the contours.

	# lower_range and upper_range by default is assigned to None. This allows
	# the contours to be altered without necessarily having to alter the 
	# lower_range and upper_range at the same time.
	def generateContours(self, contours = None, lower_range = None, 
			upper_range = None):
		"""
		Set the number of contours to generate and its range.	

		@type contours: Number
		@param contours: Number of contours to generate
		@type lower_range: Number
		@param lower_range: Lower range of contour values
		@type upper_range: Number
		@param upper_range: Upper range of contours values
		"""
	
		if(contours != None): # True if the contours is specified.
			self.__contours = contours
		if(lower_range != None): # True if the lower_range is specified.
			self.__lower_range = lower_range
		if(upper_range != None): # True if the upper_range is specified.
			self.__upper_range = upper_range

	def _generateContours(self):
		"""
		Generate the specified number of contours within the specified range.

		@attention: In order to generate an iso surface, the 'lower_range' and 
		'upper_range' must be equal.
		"""

		self.__vtk_contour.GenerateValues(self.__contours, self.__lower_range, 
				self.__upper_range)

	def _getContourModuleOutput(self):
		"""
		Return the output of the contour.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_contour.GetOutput()

	def _isContoursSet(self):
		"""
		Return whether the number of contours have been specified.

		@rtype: Boolean
		@return: True or False
		"""

		if(self.__contours != None):
			return True
		else:
			return False

	def _isLowerRangeSet(self):
		"""
		Return whether the lower range has been specified.

		@rtype: Boolean
		@return: True or False
		"""

		if(self.__lower_range != None):
			return True
		else:
			return False

	def _isUpperRangeSet(self):
		"""
		Return whether the upper range has been specified.

		@rtype: Boolean
		@return: True or False
		"""

		if(self.__upper_range != None):
			return True
		else:
			return False
