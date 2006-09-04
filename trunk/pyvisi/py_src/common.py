"""
Class that defines the common operations invoked by the components. 
"""

import vtk

class Common:
	"""
	@author: John Ngui
	@author: Lutz Gross
	"""

	def __init__(self, open_scene, data_collector):
		"""
		Initialize all the instance variables.

		@type open_scene: L{OpenScene <openscene.OpenScene>} object
		@param open_scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>} 
			object
		@param data_collector: Source of data for visualization
		"""

		self.open_scene = open_scene
		self.data_collector = data_collector
		self.vtk_mapper = None
		self.vtk_actor = None

	def setMapper(self, component):
		"""
		Set up the mapper and its input.

		@type component: String
		@param component: Component to be mapped
		"""

		self.vtk_mapper = vtk.vtkDataSetMapper()
		eval("self.vtk_mapper.SetInput(%s)" % component)

	def setActor(self):
		"""
		Set up the actor and its mapper.
		"""

		self.vtk_actor = vtk.vtkActor()
		self.vtk_actor.SetMapper(self.vtk_mapper)

	def addActor(self):
		"""
		Add the actor to the renderer.
		"""

		self.open_scene.getRenderer().AddActor(self.vtk_actor) 

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
