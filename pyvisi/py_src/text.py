"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from style import Style

class Text:
	"""
	Class that displays 2D text.
	"""

	def __init__(self, scene):
		"""
		@type scene: L{Scene <scene.Scene>}	object
		@param scene: Scene in which components are to be added to
		"""

		self.scene = scene
		self.vtk_text_mapper = vtk.vtkTextMapper()
		self.vtk_text_actor = vtk.vtkScaledTextActor()

	def setText(self, text):
		"""
		Set up the text mapper and text to be displayed.

		@type text: String
		@param text: Text to be displayed	
		"""

		self.vtk_text_mapper.SetInput(text)
		
		self.setActor()
		self.addActor()

	def setActor(self):
		"""
		Set up the 2D text actor, its mapper and its display position.
		"""

		self.vtk_text_actor.SetMapper(self.vtk_text_mapper)
		self.setPosition(50, 20) # Default text position

	def addActor(self):
		"""
		Add the 2D text actor to the renderer.
		"""

		self.scene.getRenderer().AddActor2D(self.vtk_text_actor)

	def setPosition(self, x_coor, y_coor):
		"""
		Set the display position of the text.
		@type x_coor: Number
		@param x_coor: Coordinate along the x-axis
		@type y_coor: Number
		@param y_coor: Coordinate along the y-axis
		"""

		self.vtk_text_actor.SetDisplayPosition(x_coor, y_coor)

	def setStyle(self, style):
		"""
		Set the style of the text.
		@type style: L{Style <style.Style>} object
		@param style: Style of the text
		"""

		self.vtk_text_mapper.SetTextProperty(style.getStyle())		
	
