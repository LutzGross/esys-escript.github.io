"""
@author: John NGUI
"""

import vtk
from constant import ImageFormat

class ImageReader:
	"""
	Class that defines an image reader.
	"""

	def __init__(self, format):
		"""	
		Initialise the image reader.

		@type format: String
		@param format: Format of the image 
		"""
		print "Imgae Reader"
		self.__format = format
		self.__vtk_image_reader = self.getImageReader()

	def getImageReader(self):
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

	def setFileName(self, file_name):
		"""
		@type file_name: String
		@param file_name: Image file name from which data is to be read 
		"""
		print "image reader set file name"
		self.__vtk_image_reader.SetFileName(file_name)

	def _getOutput(self):
		print "return image reader..."
		return self.__vtk_image_reader.GetOutput()

