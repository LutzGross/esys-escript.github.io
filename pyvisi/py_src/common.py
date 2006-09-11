"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk

class Common:
	"""
	Class that defines the common operations invoked by the components. 
	"""

	def __init__(self, scene, data_collector = None):
		"""
		Initialize all the instance variables.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>} 
			object
		@param data_collector: Source of data for visualization
		"""

		self.scene = scene
		self.data_collector = data_collector
		self.vtk_mapper = None
		self.vtk_actor = None

	def setMapper(self, component, lut = None):
		"""
		Set up the mapper and its input.

		@type component: String
		@param component: Component to be mapped
		@type lut: L{BlueToRed <colormap.BlueToRed>} or 
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Color lookup table to be used by the mapper
		"""

		self.vtk_mapper = vtk.vtkDataSetMapper()
		eval("self.vtk_mapper.SetInput(%s)" % component)
		
		if(lut != None):
			self.vtk_mapper.SetLookupTable(lut.getLut())	

	def setActor(self):
		"""
		Set up the actor and its mapper.
		"""

		self.vtk_actor = vtk.vtkActor()
		self.vtk_actor.SetMapper(self.vtk_mapper)

	def setTexture(self, texture):
		"""
		Set the texture of the actor.

		@type texture: vtkTexture
		@param texture: Texture map of the image
		"""
		self.vtk_actor.SetTexture(texture)

	def addActor(self):
		"""
		Add the actor to the renderer.
		"""

		self.scene.getRenderer().AddActor(self.vtk_actor) 

	def setOpacity(self, opacity):
		"""
		Set the opacity (transparency) of the actor.

		@type opacity: Number
		@param opacity: Opacity (transparency) of the actor
		"""

		self.getProperty().SetOpacity(opacity)

	def setColor(self, red, green, blue):
		self.getProperty().SetColor(red, green, blue)

	def setRepresentation(self, representation):
		"""
		Set the representation of the actor.

		@type representation: String
		@param representation: Representation type (I{i.e. Wireframe})
		"""

		eval("self.getProperty().SetRepresentationTo%s()" % representation)
	
	def getProperty(self):
		"""
		Return the property of the actor.

		@rtype: vtkProperty
		@return: VTK property	
		"""

		return self.vtk_actor.GetProperty()

		
class Component:
     pass
