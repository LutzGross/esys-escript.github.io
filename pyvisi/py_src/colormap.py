"""
Class that defines the mapping of colors.
"""

class ColorMap:
	"""
	@author: John Ngui
	@author: Lutz Gross	
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
		self.colorMap["Yelow"] = [255, 255, 0]
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
		@return Blue(B) value from the RGB
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
		

class BlueRed:
   """
   color map with spectrum from blue to red
   """
   pass

class RedBlue:
   """
   color map with spectrum from red to blue
   """
   pass
