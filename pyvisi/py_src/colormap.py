"""
@var RED: Constant representing red
@type RED: RGB color
@var GREEN: Constant representing green
@type GREEN: RGB color 
@var BLUE: Constant representing blue 
@type BLUE: RGB color
@var BLACK: Constant representing black
@type BLACK: RBG color
@var WHITE: Constant representing white
@type WHITE: RGB color
@var YELLOW: Constant representing yellow
@type YELLOW: RGB color
@var PINK: Constant represnting pink
@type PINK: RGB color
@var ORANGE: Constant representing orange
@type ORANGE: RGB color
@var PURPLE: Constant representing purple
@type PURPLE: RGB color
@var GREY: Constant representing grey
@type GREY: RGB color 
@var BROWN: Constant representing brown
@type BROWN: RGB color 

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
		self.colorMap = {}
		self.buildColorMap()

	def buildColorMap(self):
		"""
		Build a hash that defines mapping of colors.
		"""

		self.colorMap["Red"] = [255, 0, 0]
		self.colorMap["Green"]= [0, 255, 0]
		self.colorMap["Blue"] = [0, 0, 255]
		self.colorMap["Black"] = [0, 0, 0]
		self.colorMap["White"] = [255, 255, 255]
		self.colorMap["Yellow"] = [255, 255, 0]
		self.colorMap["Pink"] = [255, 20, 147]
		self.colorMap["Orange"] = [255, 69, 0]
		self.colorMap["Purple"] = [138, 43, 226]
		self.colorMap["Grey"] = [169, 169, 169]
		self.colorMap["Brown"] = [139, 69, 19]

	def getR(self):
		"""
		Return the red(R) value from the RGB.

		@rtype: Number
		@return: Red(R) value from the RGB
		"""

		return self.colorMap[self.color_name][0]

	def getG(self):
		"""
		Return the green(G) value from the RGB.

		@rtype: Number
		@return: Green(R) value from the RGB
		"""

		return self.colorMap[self.color_name][1]

	def getB(self):
		"""
		Return the blue(B) value from the RGB.

		@rtype: Number
		@return: Blue(B) value from the RGB
		"""

		return self.colorMap[self.color_name][2]

	def getColor(self):
		"""
		Return the name of the color.
	
		@rtype: String
		@return: Name of the color
		"""

		return self.color_name

# Constants
RED = ColorMap("Red")
GREEN = ColorMap("Green")
BLUE = ColorMap("Blue")
BLACK = ColorMap("Black")
WHITE = ColorMap("White")
YELLOW = ColorMap("Yellow")
PINK = ColorMap("Pink")
ORANGE = ColorMap("Orange")
PURPLE = ColorMap("Purple")
GREY = ColorMap("Grey")
BROWN = ColorMap("Brown")
	
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



