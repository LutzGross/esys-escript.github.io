"""
@author: John NGUI
"""

import vtk
from actor import Actor2D
from constant import Viewport, Color

class Text2D(Actor2D):
	"""
	Class that defines a 2D text actor.
	"""

	def __init__(self, scene, text, viewport = Viewport.SOUTH_WEST):
		"""
		Initialise the 2D text actor.

		@type text: String
		@param text: 2D text to be displayed
		"""

		self.__scene = scene
		self.__text = text
		self.__viewport = viewport
		self._vtk_actor2D = vtk.vtkTextActor()

		self.__setupText2D()
	
	def __setupText2D(self):
		"""
		Setup the 2D text.
		"""

		self.__setInput()
		# Add the 2D text to the appropriate renderer.
		self.__scene._addActor2D(self.__viewport, self._vtk_actor2D)

	def __setInput(self):
		"""
		Set the 2D text to be displayed.
		"""

		self._vtk_actor2D.SetInput(self.__text)

	def setFontSize(self, size):
		"""
		Set the 2D text size.

		@type size: Number
		@param size: Size of the 2D text
		"""

		self._vtk_actor2D.GetTextProperty().SetFontSize(size)

	def setFontToTimes(self):
		"""
		Set the 2D text font type to times new roman.
		"""

		self._vtk_actor2D.GetTextProperty().SetFontFamilyToTimes()

	def setFontToArial(self):
		"""
		Set the 2D text font type to arial.
		"""

		self._vtk_actor2D.GetTextProperty().SetFontFamilyToArial()

	def setFontToCourier(self):
		"""
		Set the 2D text front type to courier.
		"""

		self._vtk_actor2D.GetTextProperty().SetFontFamilyToCourier()

	def setJustificationToCenter(self):
		"""
		Set the 2D text to center justification.
		"""

		self._vtk_actor2D.GetTextProperty().SetJustificationToCentered()

	def setJustificationToLeft(self):
		"""
		Set the 2D text to left justification.
		"""

		self._vtk_actor2D.GetTextProperty().SetJustificationToLeft()

	def setJustificationToRight(self):
		"""
		Set the 2D text to right justification.
		"""

		self._vtk_actor2D.GetTextProperty().SetJustificationToRight()

	def boldOn(self):
		"""
		Bold the 2D text.
		"""

		self._vtk_actor2D.GetTextProperty().BoldOn()

	def shadowOn(self):
		"""
		Apply shadow onto the 2D text to ease visibility.
		"""

		self._vtk_actor2D.GetTextProperty().ShadowOn()

	def setColor(self, color):
		"""
		Set the color of the text.

		@type color: L{Color <constant.Color>} constant
		@param color: 2D text color
		"""

		self._vtk_actor2D.GetTextProperty().SetColor(color)

