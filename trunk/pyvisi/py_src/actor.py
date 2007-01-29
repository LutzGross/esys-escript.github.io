"""
@author: John NGUI
"""

import vtk

class Actor3D:
	"""
	Class that defines a 3D actor.
	"""
	
	def __init__(self, mapper):
		"""
		Initialise the 3D actor.

		@type mapper: vtkDataSetMapper
		@param mapper: Mapped data
		"""

		self.__mapper = mapper
		self.__vtk_actor3D = vtk.vtkActor()
		self.__setMapper()

	def __setMapper(self):
		"""
		Set the mapper for the 3D actor.
		"""

		self.__vtk_actor3D.SetMapper(self.__mapper)

	def _setTexture(self, texture):
		"""
		Set the texture for the 3D actor.	
		"""

		self.__vtk_actor3D.SetTexture(texture)

	def setOpacity(self, opacity):
		"""
		Set the opacity (transparency) of the 3D actor.

		@type opacity: Number (between 0 and 1)
		@param opacity: Opacity (transparency) of the 3D actor
		"""
	
		self.__vtk_actor3D.GetProperty().SetOpacity(opacity)

	def setColor(self, color):
		"""
		Set the color of the 3D actor.

		@type color: L{Color <constant.Color>} constant
		@param color: 3D actor color
		"""

		# NOTE: Must be used before actor.GetProperty().SetColor()
		# in order for the change of color to rendered objects to take effect.
		self.__mapper.ScalarVisibilityOff()

		# NOTE: Must be used after mapper.ScalarVisibilityOff()
		# in order for the change of color rendered objects to take effect.
		self.__vtk_actor3D.GetProperty().SetColor(color) 

	def setRepresentationToWireframe(self):
		"""
		Set the representation of the 3D actor to Wireframe.
		"""

		self.__vtk_actor3D.GetProperty().SetRepresentationToWireframe()
	
	def _setLineWidth(self, line_width):
		"""
		Set the line width of the 3D actor.

		@type line_width: Number
		@param line_width: 3D actor line width
		"""

		self.__vtk_actor3D.GetProperty().SetLineWidth(line_width)		


	def _getActor3D(self):
		"""
		Return the 3D actor.

		@rtype: vtkActor
		@return: 3D actor
		"""

		return self.__vtk_actor3D


############################################################################



class Actor2D:
	"""
	Class that defines a 2D actor.
	"""

	def __init__(self, mapper):
		"""
		Initialise the 2D actor.

		@type mapper: vtkImageMapper, etc
		@param mapper: Mapped data
		"""

		self.__mapper = mapper
		self._vtk_actor2D = vtk.vtkActor2D()
		self.__setMapper()

	def __setMapper(self):
		"""
		Set the mapper for the 2D actor.	
		"""

		self._vtk_actor2D.SetMapper(self.__mapper)

	def setPosition(self, position):
		"""
		Set the position of the 2D actor. Default position is the lower left
		hand corner.

		@type position: L{LocalPosition <position.LocalPosition>} object
		@param position: 2D actor position
		"""

		self._vtk_actor2D.SetPosition(position._getLocalPosition())

	def _getActor2D(self):
		"""
		Return the 2D actor.	

		@rtype: vtkActor2D
		@return 2D actor
		"""

		return self._vtk_actor2D


class TextActor(Actor2D):
	"""
	Class that defines a 2D text actor.
	"""

	def __init__(self, text):
		"""
		Initialise the 2D text actor.

		@type text: String
		@param text: 2D text to be displayed
		"""

		self.__text = text
		self._vtk_actor2D = vtk.vtkTextActor()

		self.__setInput()
		
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

