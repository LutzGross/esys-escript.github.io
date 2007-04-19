
import vtk

class ImageReslice:

	def __init__(self, object):
		self.__object = object
		self.__vtk_image_reslice = vtk.vtkImageReslice()

		self.__setupImageReslice()

	def __setupImageReslice(self):
		self.__setInput()

	def __setInput(self):
		self.__vtk_image_reslice.SetInput(self.__object)		

	def setSize(self, size):
		if(size > 1):
			size = 1 - (size - 1)
			self.__vtk_image_reslice.SetOutputSpacing(size, size, size)
		elif(size < 1):
			size = (1 - size) + 1
			self.__vtk_image_reslice.SetOutputSpacing(size, size, size)

	def _getOutput(self):
		return self.__vtk_image_reslice.GetOutput()
