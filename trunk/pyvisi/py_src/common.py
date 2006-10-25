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
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>} 
			object
		@param data_collector: Source of data for visualization
		"""

		self.scene = scene
		self.data_collector = data_collector
		self.vtk_mapper = vtk.vtkDataSetMapper() 
		self.vtk_actor = vtk.vtkActor() 

	def setMapperInput(self, component, lut = None):
		"""
		Set up the mapper.
		@type component: String
		@param component: Component to be mapped
		@type lut: L{BlueToRed <colormap.BlueToRed>} or 
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""
		
		# Convert unstructured grid data to polygonal data.	
		vtk_geometry = vtk.vtkGeometryFilter()
		vtk_geometry.SetInput(component)

		# Compute normals to ensure consistent orientation across neighbours.
		# This results in a better object being rendered. 
		#vtk_normals = vtk.vtkPolyDataNormals()
		#vtk_normals.SetInput(vtk_geometry.GetOutput())

		#self.vtk_mapper.SetInput(vtk_normals.GetOutput())
		#self.vtk_mapper.SetInput(vtk_geometry.GetOutput())
		self.vtk_mapper.SetInput(component)

	
		# Mapper uses the customized lookup table only if it is specified. 
		# Otherwise, the default one is used.
		if(lut != None):
			self.vtk_mapper.SetLookupTable(lut.getLut())	

	def setActorTexture(self, texture):
		"""
		Set the texture of the actor.
		@type texture: vtkTexture
		@param texture: Texture map of the image
		"""

		self.vtk_actor.SetTexture(texture)

	def setActorInput(self):
		"""
		Set up the actor.
		"""

		self.vtk_actor.SetMapper(self.vtk_mapper)


	def addActor(self):
		"""
		Add the actor to the renderer.
		"""

		self.scene.getRenderer().AddActor(self.vtk_actor) 

	def setActorOpacity(self, opacity):
		"""
		Set the opacity (transparency) of the actor.
		@type opacity: Number
		@param opacity: Opacity (transparency) of the actor
		"""

		self.vtk_actor.GetProperty().SetOpacity(opacity)

	def setActorColor(self, color):
		"""
		Set the color of the actor.
		@type color: RGB list
		@param color: Color of the actor
		"""
		
		self.vtk_actor.GetProperty().SetColor(color[0], color[1],
			color[2])

	def setActorRepresentation(self, representation):
		"""
		Set the representation of the actor.
		@type representation: String
		@param representation: Actor representation type (I{i.e. Wireframe})
		"""

		eval("self.vtk_actor.GetProperty().SetRepresentationTo%s()" % 
			representation)
	

class Component:
     pass
