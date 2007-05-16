"""
@author: John NGUI
"""

import vtk

class ImageReslice:
	"""
	Class that defines an image reslice used to resize static 
	(no interaction capability) images.
	"""

	def __init__(self, object):
		"""
		Initialise the image reslice.

		@type object: vtkImageData
		@param object: Image Data
		"""

		self.__object = object
		self.__vtk_image_reslice = vtk.vtkImageReslice()

		self.__setupImageReslice()

	def __setupImageReslice(self):
		"""
		Setup the image reslice.
		"""

		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the image reslice.
		"""

		self.__vtk_image_reslice.SetInput(self.__object)		

	def setSize(self, size):
		"""
		Set the size of the image, between 0 and 2. Size 1 (one) displays the 
		image in its original size (which is the default).  

		@type size: Number
		@param size: Size of the static image
		"""

		# By default, with image reslice, the larger the output spacing, the
		# smaller the image. Similarly, the smaller the output spacing, the 
		# larger the image. This behaviour is reversed so that the larger the 
		# size the image. Similarly, the smaller the size, the smaller the
		# image.
		if(size > 1):
			size = 1 - (size - 1)
			self.__vtk_image_reslice.SetOutputSpacing(size, size, size)
		elif(size < 1):
			size = (1 - size) + 1
			self.__vtk_image_reslice.SetOutputSpacing(size, size, size)

	def _getOutput(self):
		"""
		Return the output of the image reslice.

		@rtype: vtkImageData
		@return: Image data
		"""

		return self.__vtk_image_reslice.GetOutput()
