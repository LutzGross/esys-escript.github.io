"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk

class Light:
	"""
	Class that controls the light and its settings.
	"""

	def __init__(self, scene):
		"""
		@type scene: L{OpenScene <scene.Scene>} object
		@param scene: Scene in which components are to be added to	
		"""

		self.scene = scene
		self.vtk_light = None

		self.setLight()

	def setLight(self):
		"""
		Set up the light and associate it with the renderer.
		"""
		self.vtk_light = vtk.vtkLight()
		self.scene.getRenderer().AddLight(self.vtk_light)

	def setColor(self, colorMap):
		"""
		Set the color of the light.
		
		@type colorMap: L{ColorMap <colormap.ColorMap>} object
		@param colorMap: Color of the light
		"""

		self.vtk_light.SetColor(colorMap.getR(), colorMap.getG(), 
			colorMap.getB())


	def setFocalPoint(self, position):
		"""
		Set the focal point of the light.

		@type position: L{Position <geo.Position>} object
		@param position: Light focal point position
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
		@param intensity: intensity (brightness) of the light
		"""

		self.vtk_light.SetIntensity(intensity)



