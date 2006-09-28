"""
@author: John Ngui
@author: Lutz Gross
"""	

import vtk
from constants import *

class Scene:
	"""
	Class that defines a scene in which components are displayed.
	"""

	def __init__(self, renderer, x_size = 800, y_size = 600):
		"""
		@type renderer: String
		@param renderer: Type of rendering 		
		@type x_size: Number
		@param x_size: Size of the rendering window on the x-axis
		@type y_size: Number
		@param y_size: Size of the rendering window on the y-axis
		"""	

		self.renderer = renderer
		self.x_size = x_size
		self.y_size = y_size
		self.vtk_renderer = vtk.vtkRenderer() 
		self.vtk_render_window = vtk.vtkRenderWindow() 

		self.setRenderWindow()
		#self.setRenderWindowInteractor()
		
	def saveImage(self, image_name):
		"""
		Save the rendered object as an image.		
		@type image_name: String
		@param image_name: Name of the image file 
		"""

		self.vtk_render_window.Render()
		#self.vtk_render_window.OffScreenRenderingOn()
		#self.vtk_render_window.OffScreenRenderingOff()

		# Converts the output of the render window into vtkImageData.
		vtk_window_to_image = vtk.vtkWindowToImageFilter()
		vtk_window_to_image.Modified() # Update precaution
		vtk_window_to_image.Update() # Update precaution
		vtk_window_to_image.SetInput(self.vtk_render_window)
		# Force an update to the of the output image as vtk window's 
		# modification time does not get updated automatically. 
		vtk_window_to_image.Modified() # Update precaution
		vtk_window_to_image.Update() # Update precaution
		
		# Write the image to file.
		vtk_image_writer = self.getImageWriter(self.renderer)
		vtk_image_writer.SetInput(vtk_window_to_image.GetOutput())
		vtk_image_writer.SetFileName(image_name)
		vtk_image_writer.Write()


	def setRenderWindow(self):
		"""
		Set up the renderer and rendering window.
		"""

		self.vtk_render_window.AddRenderer(self.vtk_renderer)
		# Title of the render window.
		self.vtk_render_window.SetWindowName(
			"Earth Systems Science Computational Centre Rendering Tool")
		self.vtk_render_window.SetSize(self.x_size, self.y_size)
		# Default background color is white.
		self.vtk_renderer.SetBackground(
			WHITE[0], WHITE[1], WHITE[2]) 

	def setRenderWindowInteractor(self):
		""" 
		Set up the render window interactor.
		"""

		vtk_render_window_interactor = vtk.vtkRenderWindowInteractor()
		vtk_render_window_interactor.SetRenderWindow(self.vtk_render_window)
		vtk_render_window_interactor.Initialize()
		self.vtk_render_window.Render()
		vtk_render_window_interactor.Start()
		

	def render(self):
		"""
		Render the image.
		"""

		vtk_render_window_interactor = vtk.vtkRenderWindowInteractor()
		vtk_render_window_interactor.SetRenderWindow(self.vtk_render_window)
		vtk_render_window_interactor.Initialize()
		self.vtk_render_window.Render()
		vtk_render_window_interactor.Start()


	def getImageWriter(self, renderer):
		"""
		Determine the type of renderer and return the correponding image writer.
		@type renderer: String
		@param renderer: Type of renderer
		"""

		if(renderer == "vtk_jpeg"):
			return vtk.vtkJPEGWriter()
		elif(renderer == "vtk_bmp"):
			return vtk.vtkBMPWriter()
		elif(renderer == "vtk_pnm"):
			return vtk.vtkPNMWriter()
		elif(renderer == "vtk_png"):
			return vtk.vtkPNGWriter()
		elif(renderer == "vtk_tiff"):
			return vtk.vtkTIFFWriter()
		elif(renderer == "vtk_ps"):
			return vtk.vtkPostScriptWriter()

	def getRenderer(self):
		"""
		Return the renderer.

		@rtype: vtkRenderer
		@return: VTK renderer that is used to render objects on the window
		"""

		return self.vtk_renderer
