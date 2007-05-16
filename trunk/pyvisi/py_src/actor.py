"""
@author: John NGUI
"""

import vtk

class Actor3D:
	"""
	Class that defines a 3D actor.
	"""
	
	def __init__(self):
		"""
		Initialise the 3D actor.
		"""

		self.__vtk_actor3D = vtk.vtkActor()

	def _setupActor3D(self, mapper):
		"""
		Setup the 3D actor.

		@type mapper: vtkDataSetMapper
		@param mapper: Mapped data
		"""

		self.__mapper = mapper
		self.__setMapper()
		
	def __setMapper(self):
		"""
		Set the mapper of the 3D actor.
		"""

		self.__vtk_actor3D.SetMapper(self.__mapper)

	def _setTexture(self, texture):
		"""
		Set the texture of the 3D actor.	

		@type texture: vtkTexture
		@param texture: Texture of the rendered object 
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
		@param color: Color of the 3D actor 
		"""

		# NOTE: Must be used before actor.GetProperty().SetColor()
		# in order for the change of color to take effect.
		self.__mapper.ScalarVisibilityOff()

		# NOTE: Must be used after mapper.ScalarVisibilityOff()
		# in order for the change of color to take effect.
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


###############################################################################


class Actor2D:
	"""
	Class that defines a 2D actor.
	"""

	def __init__(self):
		"""
		Initialise the 2D actor.
		"""

		self._vtk_actor2D = vtk.vtkActor2D()

	def _setupActor2D(self, mapper):
		"""
		Setup the 2D actor.

		@type mapper: vtkMapper2D
		@param mapper: Mapped data
		"""

		self.__mapper = mapper
		self.__setMapper()

	def __setMapper(self):
		"""
		Set the mapper for the 2D actor.	
		"""

		self._vtk_actor2D.SetMapper(self.__mapper)

	def setPosition(self, position):
		"""
		Set the position (XY) of the 2D actor. Default position is the lower 
		left hand corner of the window / viewport.

		@type position: L{LocalPosition <position.LocalPosition>} object
		@param position: Position of the 2D actor 
		"""

		self._vtk_actor2D.SetPosition(position._getLocalPosition())

	def _getActor2D(self):
		"""
		Return the 2D actor.	

		@rtype: vtkActor2D
		@return 2D actor
		"""

		return self._vtk_actor2D

