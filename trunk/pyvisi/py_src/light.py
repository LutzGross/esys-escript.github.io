"""
@author: John NGUI
"""

import vtk
from position import GlobalPosition
from constant import Viewport

class Light:
	"""
	Class that defines a light and its settings.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST):
		"""
		Initialise the light.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to	
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		"""

		self.__scene = scene
		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__vtk_light = vtk.vtkLight()

		self.__setupLight()

	def __setupLight(self):
		"""
		Set up the light and associate it with the renderer.
		"""
		self.__scene._addLight(self.__viewport, self.__vtk_light)

	def setColor(self, color):
		"""
		Set the light color.

		@type color: L{Color <constant.Color>} constant
		@param color: Light color
		"""

		self.__vtk_light.SetColor(color)

	def setFocalPoint(self, position):
		"""
		Set the focal point of the light.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Light focal point
		"""

		self.__vtk_light.SetFocalPoint(position._getGlobalPosition())

	def setPosition(self, position):
		"""
		Set the position of the light.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Light position
		"""

		self.__vtk_light.SetPosition(position._getGlobalPosition())

	def setAngle(self, elevation = 0, azimuth = 0):
		"""
		Set the position and focal point of the light based on elevation and 
		azimuth degree.

		@type elevation: Number
		@param elevation: Degree to rotate the light to the top and bottom
		@type azimuth: Number
		@param azimuth: Degree to rotate the camera to the left and right
		"""

		# NOTE: The elevation angle of light does not seem to suffer the same
		# constraint as the elevation angle of camera where the elevation
		# angle is constraint between 90/-90.
		self.__vtk_light.SetDirectionAngle(elevation, azimuth)

	def setIntensity(self, intensity):
		"""
		Set the intensity (brightness) of the light.

		@type intensity: Number
		@param intensity: Intensity (brightness) of the light
		"""

		self.__vtk_light.SetIntensity(intensity)

