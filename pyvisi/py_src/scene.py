"""
@author: John NGUI
"""

import vtk
from constant import Renderer, Color, Viewport

class Scene:
	"""
	Class that defines a scene in which objects are to be rendered on.	
	"""
	
	def __init__(self, renderer = Renderer.ONLINE, num_viewport = 1, 
			x_size = 1152, y_size = 864):
		"""
		Initialise the scene.

		@type renderer: String
		@param renderer: Type of rendering (i.e.Renderer.ONLINE, etc) 
		@type num_viewport: Number
		@param num_viewport: Number of viewport(s) in the scene. Either 1 or 4 
		@type x_size: Number
		@param x_size: Size of the render window along the x-axis
		@type y_size: Number
		@param y_size: Size of the render window along the y-axis
		"""

		self.__renderer = renderer
		self.__num_viewport = num_viewport
		self.__x_size = x_size
		self.__y_size = y_size
	
		self.__vtk_render_window = vtk.vtkRenderWindow()
		self.__setupScene()
		
	def __setupScene(self):
		"""
		Setup the scene.
		"""

		self.__createViewport()			
		self.__addRenderer()
		self.setBackground(Color.WHITE) # Default background color is white.

		# Default title bar.
		self.setTitleBar("Earth Systems Science Computational Centre (ESSCC)")
		self.__setSize(self.__x_size, self.__y_size)
		
		if(self.__renderer == Renderer.ONLINE): # True for online rendering.
			self.__setupOnlineRendering()
		else: # True for offline rendering.
			self.__setupOfflineRendering()


	def __createViewport(self):
		"""
		Create the viewport(s) in the scene.
		"""

		# Create the renderer(s) for the viewport(s).
		self.__vtk_renderer = [] 
		for viewport in range(0, self.__num_viewport):
			self.__vtk_renderer.append(vtk.vtkRenderer())
			
		if(self.__num_viewport == 4): 
			# Renderer for the entire scene (background to the viewports).	
			self.__vtk_renderer_background = vtk.vtkRenderer()

			# Specify the positioning of the four viewports (between 0 and 1).
			self.__vtk_renderer[Viewport.SOUTH_WEST].SetViewport(
					0.0, 0.0, 0.5, 0.5)	
			self.__vtk_renderer[Viewport.NORTH_WEST].SetViewport(
					0.0, 0.5013, 0.5, 1)
			self.__vtk_renderer[Viewport.NORTH_EAST].SetViewport(
					0.501, 0.5013, 1, 1)
			self.__vtk_renderer[Viewport.SOUTH_EAST].SetViewport(
					0.501, 0.0, 1.0, 0.5)
			
	def setBackground(self, color):
		"""
		Set the background color of the scene.

		@type color: L{Color <constant.Color>} constant
		@param color: Scene background color
		"""
		
		# Color the entire scene (background) black initially. 
		# This is carried out mainly to have the borders between 
		# the viewports visibly black.
		if(self.__num_viewport == 4): 
			self.__vtk_renderer_background.SetBackground(Color.BLACK)

		for viewport in range(0, self.__num_viewport):
			self.__vtk_renderer[viewport].SetBackground(color)

	def __addRenderer(self):
		"""
		Add the renderer(s) to the render window.
		"""
		# Renderer for the black scene (background).
		if(self.__num_viewport == 4): 
			self.__vtk_render_window.AddRenderer(
					self.__vtk_renderer_background)

		for viewport in range(0, self.__num_viewport):
			self.__vtk_render_window.AddRenderer(self.__vtk_renderer[viewport])

	def setTitleBar(self, text):
		"""
		Set the text on the title bar of the render window.

		@type text: String
		@param text: Text on the title bar
		"""

		self.__vtk_render_window.SetWindowName(text)

	def __setSize(self, x_size, y_size):
		"""
		Set the size of the render window.

		@type x_size: Number
		@param x_size: Size of the render window along the x-axis
		@type y_size: Number
		@param y_size: Size of the render window along the y-axis
		"""

		self.__vtk_render_window.SetSize(x_size, y_size)	

	def __setupOnlineRendering(self):
		"""
		Setup the window interactor for online rendering.
		"""

		# Associate the window interactor with the render window.
		self.__vtk_render_window_interactor = vtk.vtkRenderWindowInteractor(					)
		self.__vtk_render_window_interactor.SetRenderWindow(
				self.__vtk_render_window)
		self.__vtk_render_window_interactor.Initialize()

	def __setupOfflineRendering(self):
		"""
		Setup the window to image filter for offline rendering.
		"""

		# Enable offscreen rendering.
		self.__vtk_render_window.OffScreenRenderingOn()
		# Convert the output of the render window into an image.
		self.__vtk_window_to_image = vtk.vtkWindowToImageFilter()
		self.__vtk_window_to_image.SetInput(self.__vtk_render_window)
		self.__vtk_image_writer = self.__getImageWriter()

	def __getImageWriter(self):
		"""
		Return the appropriate image writer based on the specified offline
		renderer.

		@rtype: vtkImageWriter
		@return: Image writer
		"""

		if(self.__renderer == Renderer.OFFLINE_JPG):
			return vtk.vtkJPEGWriter() 
		elif(self.__renderer == Renderer.OFFLINE_BMP):
			return vtk.vtkBMPWriter() 
		elif(self.__renderer == Renderer.OFFLINE_PNM):
			return vtk.vtkPNMWriter()
		elif(self.__renderer == Renderer.OFFLINE_PNG):
			return vtk.vtkPNGWriter()
		elif(self.__renderer == Renderer.OFFLINE_TIF):
			return vtk.vtkTIFFWriter()
		elif(self.__renderer == Renderer.OFFLINE_PS):
			return vtk.vtkPostScriptWriter()
	
	def saveImage(self, image_name):
		"""
		Save the rendered object as an image.

		@type image_name: String
		@param image_name: Name of the saved image.
		"""
		
		# NOTE: Render and Modified must be called everytime before writing 
		# an image. Otherwise, only the first image will always be saved.
		# This is due to the architecture of VTK.
		self.__vtk_render_window.Render()
		self.__vtk_window_to_image.Modified()
		
		# Retrieve rendered object from the window and convert it into an 
		# image.
		self.__vtk_image_writer.SetInput(
				self.__vtk_window_to_image.GetOutput())
		self.__vtk_image_writer.SetFileName(image_name)
		self.__vtk_image_writer.Write() 	

	def render(self):
		"""
		Render the object onto the scene.
		"""	

		self.__vtk_render_window.Render()

		# NOTE: Once Start is executed, the driver will not further execute 
		# any subsequent codes thereafter.
		self.__vtk_render_window_interactor.Start()
	
	def _addActor3D(self, viewport, actor):
		"""
		Add the actor3D to the appropriate viewport.

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport which the actor3D is to be added to 
		@type actor: L{Actor3D <actor.Actor3D>} object
		@param actor: Actor3D which is to be added to the viewport 
		"""

		self.__vtk_renderer[viewport].AddActor(actor)	

	def _addActor2D(self, viewport, actor):
		"""
		Add the actor2D to the appropriate viewport.

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport which the actor2D is to be added to 
		@type actor: L{Actor2D <actor.Actor2D>} object
		@param actor: Actor2D which is to be added to the viewport 
		"""

		self.__vtk_renderer[viewport].AddActor2D(actor)	

	def _setActiveCamera(self, viewport, camera):
		"""
		Set the camera to the appropriate viewport.	

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport which the camera is to be added to 
		@type camera: L{Camera <camera.Camera>} object
		@param camera: Camera which is to be assigned to the viewport
		"""

		self.__vtk_renderer[viewport].SetActiveCamera(camera)

	def _addLight(self, viewport, light):
		"""
		Add the light to the appropriate viewport.

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport which the camera is to be added to 
		@type light: L{Light <light.Light>} object
		@param light: Light which is to be assigned to the viewport
		"""

		self.__vtk_renderer[viewport].AddLight(light)

	def _getRenderer(self):
		"""
		Return the renderer(s)

		@rtype: One or more column tuple of vtkRenderer
		@return: Renderer(s)
		"""
	
		return self.__vtk_renderer
