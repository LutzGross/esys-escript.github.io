
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from constant import Renderer, Color, Viewport, ImageFormat
from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class Scene:
	"""
	Class that defines a scene. A scene is a window in which objects are to 
	be rendered on. Only one scene needs to be created. However, a scene may 
	be divided into four smaller windows called viewports (if needed). 
	Each viewport can render a different object.  

	@attention: If the IV or VRML renderer type is used, then only one 
	viewport can be used (four viewports are not supported for these two 
	cases). 
	"""

	def __init__(self, renderer = Renderer.ONLINE, num_viewport = 1, 
			x_size = 1152, y_size = 864):
		"""
		Initialise the scene.

		@type renderer: L{Renderer <constant.Renderer>} constant
		@param renderer: Type of renderer 
		@type num_viewport: Number
		@param num_viewport: Number of viewport(s) in the scene. Either 1 or 4 
		@type x_size: Number
		@param x_size: Size of the render window along the x-axis
		@type y_size: Number
		@param y_size: Size of the render window along the y-axis
		"""
                if getMPISizeWorld()>1:
                   raise ValueError,"pyvisi.Scene is not running on more than one processor."

		self.__renderer = renderer
		self.__num_viewport = num_viewport
		self.__x_size = x_size
		self.__y_size = y_size
		# Stores the visualization modules (i.e. Map, Velocity, Ellipsoid, etc)
		# which needs to be rendered.
		self.__visualization_modules = []

		self.__OFFLINE = "offline"
		self.__VRML = "vrml"
		self.__IV = "iv"
	
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

		# True for Online rendering.
		if(self.__renderer.startswith(Renderer.ONLINE)): 
			self.__setupOnlineRendering()
			# True for all Online renderers except Renderer.ONLINE.
			if(self.__renderer != Renderer.ONLINE):
				self.__setupWindowToImage()
		# True for Offline rendering.
		elif(self.__renderer.startswith(self.__OFFLINE)): 
			self.__setupOfflineRendering()
			self.__setupWindowToImage()
		# True for Display rendering.
		elif(self.__renderer.startswith(Renderer.DISPLAY)): 
			# True for all Display renderers except Renderer.DISPLAY.
			if(self.__renderer != Renderer.DISPLAY):
				self.__setupWindowToImage()

	def __createViewport(self):
		"""
		Create the viewport(s) in the scene.
		"""

		# Create the renderer(s) for the viewport(s).
		self.__vtk_renderer = [] 
		for viewport in range(0, self.__num_viewport):
			self.__vtk_renderer.append(vtk.vtkRenderer())
			
		if(self.__num_viewport == 4): 
			# Renderer for the entire scene (background).	
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

		# Add the renderer for the black scene (background).
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
		self.__vtk_render_window_interactor = vtk.vtkRenderWindowInteractor()
		self.__vtk_render_window_interactor.SetRenderWindow(
				self.__vtk_render_window)
		self.__vtk_render_window_interactor.Initialize()

	def __setupOfflineRendering(self):
		"""
		Enables the offline rendering (no window comes up).
		"""

		# Enables the offscreen rendering.
		self.__vtk_render_window.OffScreenRenderingOn()

	def __setupWindowToImage(self):
		"""	
		Setup the window to image filter to convert the output from the render 
		window into an image, vrml or open inventor file.
		"""

		self.__vtk_window_to_image = vtk.vtkWindowToImageFilter()
		self.__vtk_window_to_image.SetInput(self.__vtk_render_window)
		self.__vtk_image_writer = self.__getImageWriter()
		
	def __getImageWriter(self):
		"""
		Return the appropriate image writer or exporter based on the 
		specified renderer.

		@rtype: vtkImageWriter or vtkExporter
		@return: Image writer or exporter
		"""

		if(self.__renderer.endswith(ImageFormat.JPG)):
			return vtk.vtkJPEGWriter() 
		elif(self.__renderer.endswith(ImageFormat.BMP)):
			return vtk.vtkBMPWriter() 
		elif(self.__renderer.endswith(ImageFormat.PNM)):
			return vtk.vtkPNMWriter()
		elif(self.__renderer.endswith(ImageFormat.PNG)):
			return vtk.vtkPNGWriter()
		elif(self.__renderer.endswith(ImageFormat.TIF)):
			return vtk.vtkTIFFWriter()
		elif(self.__renderer.endswith(ImageFormat.PS)):
			return vtk.vtkPostScriptWriter()
		elif(self.__renderer.endswith(self.__VRML)):
			return vtk.vtkVRMLExporter() # Generates VRML files (.wrl).
		elif(self.__renderer.endswith(self.__IV)):
			return vtk.vtkIVExporter() # Generate OpenInventor files (.iv).
	
	def __saveImage(self, image_name):
		"""
		Save the rendered object as an image, vrml or open inventor file.

		@type image_name: String
		@param image_name: Name of the saved image, vrml or open inventor file
		"""

		# NOTE: Render and Modified must be called everytime before writing 
		# an image, vrml or open inventor file. Otherwise, only the first 
		# object will always be saved. This is due to the architecture of VTK.
		self.__vtk_render_window.Render()
		self.__vtk_window_to_image.Modified()
		
		# Retrieve the rendered object from the window and convert it into an 
		# image, vrml or open inventor file.
		# True for all writers besides VRML.
		if(not(self.__renderer.endswith(self.__VRML)) and \
				not(self.__renderer.endswith(self.__IV))): 
			self.__vtk_image_writer.SetInput(
					self.__vtk_window_to_image.GetOutput())
		# True only for VRML and IV.
		elif(self.__renderer.endswith(self.__VRML) or \
				self.__renderer.endswith(self.__IV)):
			self.__vtk_image_writer.SetInput(
					self.__vtk_render_window)
		self.__vtk_image_writer.SetFileName(image_name)
		self.__vtk_image_writer.Write() 	

	def __animate(self):
		"""	
		Animate the rendered object on-the-fly.
		"""

		# With Render() ONLY, the rendered object is animated onto the 
		# scene on-the-fly and no interaction can occur.
		self.__vtk_render_window.Render()

	def render(self, image_name = None):
		"""
		Render the object using either the online, offline or display mode.

		@type image_name: String
		@param image_name: Name of the saved image, vrml or open inventor file
		"""	

		for i in range(0, len(self.__visualization_modules)):
			self.__visualization_modules[i]._render(self)	
			self.__vtk_render_window.Render()

		self.__vtk_render_window.Render()

		if(self.__renderer.startswith(Renderer.ONLINE)):
			# NOTE: Once Start() is executed, the driver will not further 
			# execute any subsequent codes thereafter unless the 'q' or 
			# 'e' keys are pressed.
			self.__vtk_render_window_interactor.Start()

			# True for all online renderers except Renderer.ONLINE.
			if(self.__renderer != Renderer.ONLINE):
				self.__saveImage(image_name)					
		# True for all display renderers except Renderer.DISPLAY.
		elif(self.__renderer.startswith(self.__OFFLINE) or 
				self.__renderer != Renderer.DISPLAY):
			self.__saveImage(image_name)					
	
	def _addActor3D(self, viewport, actor):
		"""
		Add the actor3D to the appropriate viewport.

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport in which the actor3D is to be added to 
		@type actor: vtkActor
		@param actor: Actor3D which is to be added to the viewport 
		"""

		self.__vtk_renderer[viewport].AddActor(actor)	

	def _addActor2D(self, viewport, actor):
		"""
		Add the actor2D to the appropriate viewport.

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport in which the actor2D is to be added to 
		@type actor: vtkActor2D
		@param actor: Actor2D which is to be added to the viewport 
		"""

		self.__vtk_renderer[viewport].AddActor2D(actor)	

	def _setActiveCamera(self, viewport, camera):
		"""
		Set the camera to the appropriate viewport.	

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport in which the camera is to be added to 
		@type camera: vtkCamera
		@param camera: Camera which is to be assigned to the viewport
		"""

		self.__vtk_renderer[viewport].SetActiveCamera(camera)

	def _addLight(self, viewport, light):
		"""
		Add the light to the appropriate viewport.

		@type viewport: L{Viewport <constant.Viewport>} constant 
		@param viewport: Viewport in which the camera is to be added to 
		@type light: vtkLight
		@param light: Light which is to be assigned to the viewport
		"""

		self.__vtk_renderer[viewport].AddLight(light)

	def _getRenderer(self):
		"""
		Return the renderer(s)

		@rtype: List 
		@return: A list of renderer(s)
		"""
	
		return self.__vtk_renderer

	def _addVisualizationModules(self, module):
		"""
		Store visualization modules (i.e. Map, Velocity, Ellipsoid, etc) 
		which needs to be rendered.
		
		@type module: Visualization modules
		@param module: Visualization modules to be rendered
		"""

		self.__visualization_modules.append(module)

