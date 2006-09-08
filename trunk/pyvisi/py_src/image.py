import vtk
from common import *

class Image(Common):

	def __init__(self, scene, format):
		Common.__init__(self, scene)
		self.vtk_image_reader = self.determineReader(format) 
		self.vtk_texture = None
		self.vtk_plane = None

	def determineReader(self, format):
		if(format == "jpg"):
			return vtk.vtkJPEGReader()	
		elif(format == "bmp"):
			return vtk.vtkBMPReader()

	def setFileName(self, file_name):
		self.vtk_image_reader.SetFileName(file_name)

		self.setTexture()
		self.setPlane()

		Common.setMapper(self, "self.vtk_plane.GetOutput()")
		Common.setActor(self)
		Common.setTexture(self.vtk_texture)
		Common.addActor(self)

	def setTexture(self):
		self.vtk_texture = vtk.vtkTexture()
		self.vtk_texture.SetInput(self.vtk_image_reader.GetOutput())	

	def setPlane(self):
		self.vtk_plane = vtk.vtkPlaneSource()
		
		
		
		
		
		


