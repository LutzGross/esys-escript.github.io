"""
@author: John NGUI
"""

import vtk
from position import GlobalPosition
from constant import Viewport

class Camera:
	"""
	Class that defines a camera and its settings.	
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
		@param data_collector: Deal with source of data for visualisation
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
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

		center = self.__data_collector._getOutput().GetCenter()
		# Default camera focal point is the center of the object.
		self.setFocalPoint(GlobalPosition(center[0], center[1], center[2]))
		# Default camera position is the center of the object but with a slight
		# distance on the z-axis.
		self.setPosition(GlobalPosition(center[0], center[1], center[2] * 3))
		# Assign the camera to the appropriate renderer
		#self.__scene._getRenderer()[self.__viewport].SetActiveCamera(
		#		self.__vtk_camera)
		self.__scene._setActiveCamera(self.__viewport, self.__vtk_camera)
		self.__resetCamera()	

	def setFocalPoint(self, position):
		"""
		Set the focal point of the camera.
		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera focal point
		"""

		self.__vtk_camera.SetFocalPoint(position._getGlobalPosition())		
		self.__resetCamera()

	def setPosition(self, position):
		"""
		Set the position of the camera.
		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera position
		"""

		self.__vtk_camera.SetPosition(position._getGlobalPosition())
		self.__resetCamera()

	def setClippingRange(self, near_clipping, far_clipping):
		"""
		Set the near and far clipping plane of the camera.
		@type near_clipping: Number
		@param near_clipping: Distance to the near clipping plane
		@type far_clipping: Number
		@param far_clipping: Distance to the far clipping plane
		"""
	
		self.vtk__camera.SetClippingRange(near_clipping, far_clipping)
		self.__resetCamera()

	def setViewUp(self, position):
		"""
		Set the view up direction of the camera.
		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera view up direction
		"""

		self.__vtk_camera.SetViewUp(position._getGlobalPosition())
		self.__resetCamera()

	def azimuth(self, angle):
		"""
		Rotate the camera to the left and right.
		@type angle: Number
		@param angle: Degree to rotate the camera
		"""

		self.__vtk_camera.Azimuth(angle)
		self.__resetCamera()		

	def elevation(self, angle):
		"""
		Rotate the camera to the top and bottom.
		@type angle: Number
		@param angle: Degree to rotate the camera (only between -90 and 90)
		"""

		self.__vtk_camera.Elevation(angle)
		# Recompute the view up vector. If not used the elevation angle is  
		# unable to exceed 87/-87 degrees. Secondly, a warning resetting the
		# view up will also be thrown and the rendered object may be incorrect.
		# With the view up recomputed, the elevation angle can reach between
		# 90/-90 degress. Exceeding that, the rendered object may be incorrect.
		self.__vtk_camera.OrthogonalizeViewUp()
		self.__resetCamera()

	def roll(self, angle):
		"""
		Roll the camera to the left and right.
		@type angle: Number
		@param angle: Degree to turn the camera
		"""

		self.__vtk_camera.Roll(-angle)

	def backView(self):
		"""
		View the back of the rendered object.
		"""

		self.azimuth(180)
		self.__resetCamera()

	def topView(self):
		"""
		View the top of the rendered object.
		"""
		
		self.elevation(90)
		self.__resetCamera()

	def bottomView(self):
		"""
		View the bottom of the rendered object.
		"""

		self.elevation(-90)
		self.__resetCamera()

	def leftView(self):
		"""
		View the left side of the rendered object.
		"""

		self.azimuth(-90)
		self.__resetCamera()

	def rightView(self):
		"""
		View the right side of the rendered object.
		"""

		self.azimuth(90)
		self.__resetCamera()

	def isometricView(self):
		"""
		View the isometric side of the rendered object.
		"""

		self.roll(-45)
		self.elevation(-45)
		
	def __resetCamera(self):
		"""
		Reposition the camera so that all actors can be seen. Needs to
		be called whenever the camera's settings are modified.
		"""

		self.__scene._getRenderer()[self.__viewport].ResetCamera() 
		
