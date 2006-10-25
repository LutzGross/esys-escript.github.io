"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common 

class Image(Common):
	"""
	Class that displays an image.
	"""

	def __init__(self, scene, format):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type format: String
		@param format: Format of the image (i.e. jpeg)
		"""

		Common.__init__(self, scene)
		self.vtk_image_reader = self.getImageReader(format) 
		self.vtk_texture = vtk.vtkTexture()

	def getImageReader(self, format):
		"""
		Determine the image format and return the corresponding image reader.
		@type format: String
		@param format: Format of the image 
		@rtype: vtkImageReader2 (i.e. vtkJPEGReader)
		@return: VTK image reader that is used to read an image
		"""

		if(format == "jpeg"):
			return vtk.vtkJPEGReader()	
		elif(format == "bmp"):
			return vtk.vtkBMPReader()
		elif(format == "pnm"):
			return vtk.vtkPNMReader()
		elif(format == "png"):
			return vtk.vtkPNGReader()
		elif(format == "tiff"):
			return vtk.vtkTIFFReader()

	def setFileName(self, file_name):
		"""
		Set the file name, setup the mapper and the actor.	
		@type file_name: String
		@param file_name: Image file name from which data is to be read 
		"""

		vtk_plane = vtk.vtkPlaneSource()
		self.vtk_image_reader.SetFileName(file_name)
		self.setTexture()

		Common.setMapperInput(self, vtk_plane.GetOutput())
		Common.setActorTexture(self, self.vtk_texture)
		Common.setActorInput(self)
		Common.addActor(self)

	def setTexture(self):
		"""
		Set the texture map.	
		"""

		self.vtk_texture.SetInput(self.vtk_image_reader.GetOutput())	


