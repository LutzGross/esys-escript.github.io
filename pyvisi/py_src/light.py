"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk

class Light:
	"""
	Class that controls the light and its settings.
	"""

	def __init__(self, scene, data_collector):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to	
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		"""

		self.scene = scene
		self.data_collector = data_collector
		self.vtk_light = vtk.vtkLight()

		self.setLight()

	def setLight(self):
		"""
		Set up the light and associate it with the renderer.
		"""
		self.scene.getRenderer().AddLight(self.vtk_light)

	def setColor(self, color):
		"""
		Set the light color.
		@type color: RGB list 
		@param color: Color of the light
		"""

		self.vtk_light.SetColor(color[0], color[1], color[2])

	def setFocalPoint(self, position):
		"""
		Set the focal point of the light.
		@type position: L{Position <geo.Position>} object
		@param position: Light focal point 
		"""

		self.vtk_light.SetFocalPoint(position.getXCoor(), position.getYCoor(),
			position.getZCoor())

	def setPosition(self, position):
		"""
		Set the position of the light.
		@type position: L{Position <geo.Position>} object
		@param position: Light position
		"""

		self.vtk_light.SetPosition(position.getXCoor(), position.getYCoor(),
			position.getZCoor())

	def setIntensity(self, intensity):
		"""
		Set the intensity (brightness) of the light.
		@type intensity: Number
		@param intensity: Intensity (brightness) of the light
		"""

		self.vtk_light.SetIntensity(intensity)

