"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from constants import *

class Style:
	"""
	Class that define the style of text.
	"""

	def __init__(self):
		self.vtk_text_property = vtk.vtkTextProperty() 

	def setFontFamily(self, family):
		"""
		Set the font family (i.e. Times, Arial).
		@type family: String
		@param family: Type of font
		"""
		eval("self.vtk_text_property.SetFontFamilyTo%s()" % family)

	def boldOn(self):
		"""
		Bold the text.
		"""

		self.vtk_text_property.BoldOn()

	def italicOn(self):
		"""
		Italize the text.
		"""

		self.vtk_text_property.ItalicOn()

	def shadowOn(self):
		"""
		Apply shadows on the text.
		"""

		self.vtk_text_property.ShadowOn()

	def setColor(self, color):
		"""
		Set the text color.
		@type color: RGB list
		@param color: Color of the text
		"""

		self.vtk_text_property.SetColor(color[0], color[1],
			color[2])

	def getStyle(self):
		"""
		Return the style object.
		@rtype: vtkTextProperty
		@return: VTK text property that is used to specify text style
		"""

		return self.vtk_text_property

