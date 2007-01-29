"""
@author: John NGUI
"""

import vtk

class Texture:
	"""
	Class that defines a texture map.
	"""

	def __init__(self, image):
		"""
		Initialise the texture.

		@type: vtkImageData
		@param: Image data from which data is read
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
		@return: Texture of the map
		"""

		return self.__vtk_texture

