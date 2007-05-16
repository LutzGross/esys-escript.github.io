"""
@author: John NGUI
"""

import vtk
from constant import ImageFormat

class ImageReader:
	"""
	Class that defines an image reader. An image reader is used to read
	data from an image in a variety of formats.
	"""

	def __init__(self, format):
		"""	
		Initialise the image reader.

		@type format:  L{ImageFormat <constant.ImageFormat>} constant
		@param format: Format of the image 
		"""

		self.__format = format
		self.__vtk_image_reader = self.__getImageReader()

	def __getImageReader(self):
		"""
		Return the corresponding image reader based on the supplied image 
		format.

		@rtype: vtkImageReader2 (i.e. vtkJPEGReader, etc)
		@return: Image reader 
		"""

		if(self.__format == ImageFormat.JPG):
			return vtk.vtkJPEGReader()	
		elif(self.__format == ImageFormat.BMP):
			return vtk.vtkBMPReader()
		elif(self.__format == ImageFormat.PNM):
			return vtk.vtkPNMReader()
		elif(self.__format == ImageFormat.PNG):
			return vtk.vtkPNGReader()
		elif(self.__format == ImageFormat.TIF):
			return vtk.vtkTIFFReader()

	def setImageName(self, image_name):
		"""
		Set the image file name to be read.

		@type image_name: String
		@param image_name: Image name from which data is to be read 
		"""

		self.__vtk_image_reader.SetFileName(image_name)

	def _getImageReaderOutput(self):
		"""
		Return the output of the image reader.

		@rtype: vtkImageData
		@return: Image data 
		"""

		return self.__vtk_image_reader.GetOutput()

