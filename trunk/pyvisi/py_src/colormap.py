"""
@author: John Ngui
@author: Lutz Gross	
"""

class ColorMap:
	"""
	Class that defines the mapping of colors.
	"""

	def __init__(self, color_name):
		"""
		@type color_name: String
		@param color_name: Name of the color
		"""

		self.color_name = color_name
		self.color_map = {}
		self.buildColorMap()

	def buildColorMap(self):
		"""
		Build a hash that defines mapping of colors.
		"""

		self.color_map["Red"] = [1, 0, 0]
		self.color_map["Green"]= [0, 1, 0]
		self.color_map["Blue"] = [0, 0, 1]
		self.color_map["Black"] = [0, 0, 0]
		self.color_map["White"] = [1, 1, 1]
		self.color_map["Yellow"] = [1,1,0]
		self.color_map["Pink"] = [1, 0.0784, 0.4588]
		self.color_map["Orange"] = [1, 0.2706, 0]
		self.color_map["Purple"] = [0.5412, 0.1680, 0.8828]
		self.color_map["Grey"] = [0.6602, 0.6602, 0.6602]
		self.color_map["Brown"] = [0.5430, 0.2700, 0.0742]

	def getR(self):
		"""
		Return the red(R) value from the RGB.

		@rtype: Number
		@return: Red(R) value from the RGB
		"""

		return self.color_map[self.color_name][0]

	def getG(self):
		"""
		Return the green(G) value from the RGB.

		@rtype: Number
		@return: Green(G) value from the RGB
		"""

		return self.color_map[self.color_name][1]

	def getB(self):
		"""
		Return the blue(B) value from the RGB.

		@rtype: Number
		@return: Blue(B) value from the RGB
		"""

		return self.color_map[self.color_name][2]

	def getColor(self):
		"""
		Return the name of the color.
	
		@rtype: String
		@return: Name of the color
		"""

		return self.color_name
	
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



