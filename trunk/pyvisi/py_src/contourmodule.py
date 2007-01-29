"""
@author: John NGUI
"""

import vtk

class ContourModule:
	"""
	Class that defines contour.
	"""

	def __init__(self, object):
		"""
		Initliase the contour.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the contour
		"""

		self.__object = object
		self.__vtk_contour = vtk.vtkContourFilter()

		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the contour.
		"""

		self.__vtk_contour.SetInput(self.__object)

	# lower_range and upper_range by default is assigned to None. This allows
	# the contours to be altered without necessarily having to alter the 
	# lower_range and upper_range.
	def generateContours(self, contours, lower_range = None, 
			upper_range = None):
		"""
		Generate the specified number of contours withing the specified range.

		@type contours: Number
		@param contours: Number of contours to generate
		@type lower_range: Number
		@param lower_range: Lower range of the contours
		@type upper_range: Number
		@param upper_range: Upper range of the contours
		"""

		if(lower_range != None): # True if the lower_range is altered.
			self.__lower_range = lower_range
		if(upper_range != None): # True if the upper_range is altered.
			self.__upper_range = upper_range

		self.__vtk_contour.GenerateValues(contours, self.__lower_range, 
				self.__upper_range)


	# NOTE: Method not used at the moment.
	def setValue(self, contour, value):
		"""
		Set the specific contour value.

		@type contour: Number
		@param contour: Contour 
		@type value: Number
		@param value: Contour value
		"""

		self.__vtk_contour.SetValue(contour, value)
		

	def _getOutput(self):
		"""
		Return the output of the contour.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_contour.GetOutput()

