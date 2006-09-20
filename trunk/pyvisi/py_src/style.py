import vtk
from constants import *

class Style:

	def __init__(self):
		self.vtk_text_property = vtk.vtkTextProperty() 

	def setFontFamily(self, family):
		eval("self.vtk_text_property.SetFontFamilyTo%s()" % family)

	def setBold(self):
		self.vtk_text_property.BoldOn()

	def setItalic(self):
		self.vtk_text_property.ItalicOn()

	def setShadow(self):
		self.vtk_text_property.ShadowOn()

	def setColor(self, color):
		self.vtk_text_property.SetColor(color[0], color[1],
			color[2])

	def getStyle(self):
		return self.vtk_text_property

