"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


import vtk

class LookupTable:
	"""
	Class that defines a lookup table for mapping scalar values into colors.
	"""

	def __init__(self):
		"""
		Initialise the lookup table.
		"""

		self.__vtk_lookup_table = vtk.vtkLookupTable()
		self.__vtk_inverse_lookup_table = vtk.vtkLookupTable()
		self.__build()

	def __build(self):
		"""
		Generates the lookup table.
		"""

		# NOTE: Build have to be executed prior to using SetTableValue (if any).
		self.__vtk_lookup_table.Build()
		self.__vtk_inverse_lookup_table.Build()

	def _setTableValue(self):
		"""
		Setup the lookup table with colors.
		"""

		# NOTE: The color values are inversed because VTK's default lookup
		# table is inversed by itself. SetTableValue have to be executed after
		# the Build.
		for i in range(256):
			self.__vtk_lookup_table.SetTableValue(
					i, self.__vtk_inverse_lookup_table.GetTableValue(255 - i))

	def _setLookupTableToGreyScale(self):
		"""
		Setup the lookup table with grey scale.	
		"""

		self.__setHueRange(0,0)
		self.__setSaturationRange(0,0)
		self.__setValueRange(1,0)
		self.__setNumberOfTableValues(256)
		self.__setRampToSQRT()

	def __setValueRange(self, lower_range, upper_range):
		"""
		Set the value range (brightness) for the lookup table (between 0 and 1).

		@type lower_range: Number
		@param lower_range:Lower value range 
		@type upper_range: Number
		@param upper_range: Upper value range
		"""

		self.__vtk_lookup_table.SetValueRange(lower_range, upper_range)

	def __setHueRange(self, lower_range, upper_range):
		"""
		Set the hue (color) range for the lookup table (between 0 and 1).

		@type lower_range: Number
		@param lower_range:Lower hue range 
		@type upper_range: Number
		@param upper_range: Upper hue range
		"""

		self.__vtk_lookup_table.SetHueRange(lower_range, upper_range)

	def __setSaturationRange(self, lower_range, upper_range):
		"""
		Set the saturation (vibrancy) range for the lookup table \
		(between 0 and 1).

		@type lower_range: Number
		@param lower_range:Lower saturantion range 
		@type upper_range: Number
		@param upper_range: Upper saturation range
		"""

		self.__vtk_lookup_table.SetSaturationRange(lower_range, upper_range)
	
	def __setRampToSQRT(self):
		"""
		Set the table ramp to SQRT. The default ramp is S-curve.	
		"""

		self.__vtk_lookup_table.SetRampToSQRT()

	def __setNumberOfTableValues(self, table_values):
		"""
		Set the number of values (i.e. colors) in the lookup table.	

		@type table_values: Number
		@param table_values: Number of colors in the lookup table.
		"""

		self.__vtk_lookup_table.SetNumberOfTableValues(table_values)

	def _getLookupTable(self):
		"""
		Return the lookup table.

		@rtype: vtkLookupTable
		@return: Lookup table
		"""

		return self.__vtk_lookup_table
	


