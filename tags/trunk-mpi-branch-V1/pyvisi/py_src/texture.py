"""
@author: John NGUI
"""

import vtk

class Texture:
	"""
	Class that defines the texture for the rendered object.
	"""

	def __init__(self, image):
		"""
		Initialise the texture.

		@type image: vtkImageData
		@param image: Image from which data is to be read
		"""

		self.__image = image
		self.__vtk_texture = vtk.vtkTexture()

		self.__setInput()	

	def __setInput(self):
		"""
		Set the input for the texture.
		"""

		self.__vtk_texture.SetInput(self.__image)

	def _getTexture(self):
		"""
		Return the texture.

		@rtype: vtkTexture
		@return: Texture of the rendered object
		"""

		return self.__vtk_texture

