"""
@author: John NGUI
"""

import vtk
from position import GlobalPosition
from constant import Viewport

class Camera:
	"""
	Class that defines a camera. A camera controls the display angle of
	the rendered object and one is usually created for a
	L{Scene <scene.Scene>}. However, if a L{Scene <scene.Scene>} has four
	viewports, then a separate camera may be created for each viewport.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST):
		"""
		Initialise the camera.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with the source data for vizualisation
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		"""

		self.__scene = scene
		self.__data_collector = data_collector
		self.__viewport = viewport

		self.__vtk_camera = vtk.vtkCamera()
		self.__setupCamera()
		
	def __setupCamera(self):
		"""
		Setup the camera.
		"""

		# Assign the camera to the appropriate renderer
		self.__scene._setActiveCamera(self.__viewport, self.__vtk_camera)
		self.__resetCamera()	

	def setFocalPoint(self, position):
		"""
		Set the focal point of the camera.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera focal point
		"""

		self.__vtk_camera.SetFocalPoint(position._getGlobalPosition())		
		self.__resetCameraClippingRange()

	def setPosition(self, position):
		"""
		Set the position of the camera.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera position
		"""

		self.__vtk_camera.SetPosition(position._getGlobalPosition())
		self.__resetCameraClippingRange()

	def setClippingRange(self, near_clipping, far_clipping):
		"""
		Set the near and far clipping plane of the camera.

		@type near_clipping: Number
		@param near_clipping: Distance to the near clipping plane
		@type far_clipping: Number
		@param far_clipping: Distance to the far clipping plane
		"""
	
		self.__vtk_camera.SetClippingRange(near_clipping, far_clipping)

	def setViewUp(self, position):
		"""
		Set the view up direction of the camera.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera view up direction
		"""

		self.__vtk_camera.SetViewUp(position._getGlobalPosition())

	def azimuth(self, angle):
		"""
		Rotate the camera to the left and right.

		@type angle: Number
		@param angle: Degree to rotate the camera
		"""

		self.__vtk_camera.Azimuth(angle)
		self.__resetCameraClippingRange()

	def elevation(self, angle):
		"""
		Rotate the camera to the top and bottom.
		
		@type angle: Number
		@param angle: Degree to rotate the camera (only between -90 and 90)
		"""

		self.__vtk_camera.Elevation(angle)
		# Recompute the view up vector. If not used the elevation angle is  
		# unable to exceed 87/-87 degrees. A warning resetting the
		# view up will also be thrown and the rendered object may be incorrect.
		# With the view up recomputed, the elevation angle can reach between
		# 90/-90 degrees. Exceeding that, the rendered object may be incorrect.
		self.__vtk_camera.OrthogonalizeViewUp()
		self.__resetCameraClippingRange()

	def roll(self, angle):
		"""
		Roll the camera to the left and right.

		@type angle: Number
		@param angle: Degree to roll the camera
		"""

		self.__vtk_camera.Roll(-angle)
		self.__resetCameraClippingRange()

	def backView(self):
		"""
		Rotate the camera to view the back of the rendered object.
		"""

		self.azimuth(180)

	def topView(self):
		"""
		Rotate the camera to view the top of the rendered object.
		"""
		
		self.elevation(90)

	def bottomView(self):
		"""
		Rotate the camera to view the bottom of the rendered object.
		"""

		self.elevation(-90)

	def leftView(self):
		"""
		Rotate the camera to view the left side of the rendered object.
		"""

		self.azimuth(-90)

	def rightView(self):
		"""
		Rotate the camera to view the right side of the rendered object.
		"""

		self.azimuth(90)

	def isometricView(self):
		"""
		Rotate the camera to view the isometric angle of the rendered object.
		"""

		self.roll(-45)
		self.elevation(-45)
		
	def dolly(self, distance):
		"""
		Move the camera towards (greater than 1) and away (less than 1) from 
		the rendered object. 

		@type distance: Number
		@param distance: Amount to move towards or away the rendered object
		"""

		self.__vtk_camera.Dolly(distance)
		self.__resetCameraClippingRange()

	def __resetCameraClippingRange(self):
		"""
		Reset the camera clipping range based on the bounds of the visible 
		actors. This ensures the rendered object is not cut-off.
		Needs to be called whenever the camera's settings are modified.
		"""

		self.__scene._getRenderer()[self.__viewport].ResetCameraClippingRange() 

	def __resetCamera(self):
		"""
		Repositions the camera to view the center point of the actors.
		"""

		self.__scene._getRenderer()[self.__viewport].ResetCamera() 
		
