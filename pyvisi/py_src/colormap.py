"""
@author: John Ngui
@author: Lutz Gross	
"""

import vtk	

class Lut:
	"""
	Class that provides the functions to create a map spectrum.
	"""

	def __init__(self):
		self.vtk_lut = vtk.vtkLookupTable() 
		
	def setHue(self, lower_range, upper_range):
		"""
		Set the upper and lower hue(color) range.

		@type lower_range: Number
		@param lower_range: Lower range of the hue 
		@type upper_range: Number
		@param upper_range: Upper range of the hue 
		"""

		self.vtk_lut.SetHueRange(lower_range, upper_range)

	def setSaturation(self, lower_range, upper_range):
		"""
		Set the upper and lower saturation(vibrancy) range.
	
		@type lower_range: Number
		@param lower_range: Lower range of the saturation
		@type upper_range: Number
		@param upper_range: Higher range of the saturation
		"""

		self.vtk_lut.SetSaturationRange(lower_range, upper_range)

	def setValue(self, lower_range, upper_range):
		"""
		Set the upper and lower value(brightness) range.

		@type lower_range: Number
		@param lower_range: Lower range of the value 
		@type upper_range: Number
		@param upper_range: Upper range of the value
		"""
	
		self.vtk_lut.SetValueRange(lower_range, upper_range)
		
	def getLut(self):
		"""
		Return the VTK lookup table.

		@rtype: vtkLookupTable
		@return: VTK Lookup table
		"""

		return self.vtk_lut

class BlueToRed(Lut):
	"""
	Class that creates a map with spectrum from blue to red.
	"""

	def __init__(self):
		Lut.__init__(self)
		self.setHue(0.667, 0.0)	

class RedToBlue(Lut):
	"""
	Class that creates a map with spectrum from red to blue.
	"""

	def __init__(self):
		Lut.__init__(self)
		self.setHue(0.0, 0.667)



