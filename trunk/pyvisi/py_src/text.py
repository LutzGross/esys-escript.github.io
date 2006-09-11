"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk

class Text:
	"""
	Class that displays text.
	"""

	def __init__(self, scene):
		"""
		@type scene: L{Scene <scene.Scene>}	object
		@param scene: Scene in which components are to be added to
		"""

		self.scene = scene
		self.vtk_text_mapper = None
		self.vtk_text_actor = None

	def setText(self, text):
		"""
		Setup the text mapper and its input together with the default settings.

		@type text: String
		@param text: Text to be displayed	
		"""

		self.vtk_text_mapper = vtk.vtkTextMapper()
		self.vtk_text_mapper.SetInput(text)
		
		#self.setFontSize(10)
		self.setFontFamily("Times")
		#self.setJustification("Right")
		self.bold()
		self.italic()
		self.shadow()
		# Default text color is black
		self.getProperty().SetColor(0, 0, 0)
	
		self.setActor()
		self.addActor()

	def getProperty(self):
		"""
		Return the property of the text mapper.

		@rtype: vtkTextProperty
		@return: VTK text property
		"""

		return self.vtk_text_mapper.GetTextProperty()

	#def setFontSize(self, size):
	#	eval("self.getProperty().SetFontSize(%s)" % size)

	def setFontFamily(self, font):
		"""
		Set the font of the text.

		@type font: String
		@param font: Font of the text (i.e. Times, Arial and Courier)
		"""

		eval("self.getProperty().SetFontFamilyTo%s()" % font)
	
	#def setJustification(self, justification):
	#	eval("self.getProperty().SetJustificationTo%s()" % justification) 

	def bold(self):
		"""
		Bold the text.
		"""

		self.getProperty().BoldOn()

	def italic(self):
		"""
		Italize the text.
		"""

		self.getProperty().ItalicOn()

	def shadow(self):
		"""
		Shadow the text.
		"""

		self.getProperty().ShadowOn()

	def setTextColor(self, colorMap):
		"""
		Set the color of the text.

		@type colorMap: L{ColorMap <colormap.ColorMap>} object
		@param colorMap: Color of the text
		"""

		#self.vtk_text_actor.GetProperty().SetColor(colorMap.getR(),
		#	colorMap.getG(), colorMap.getB())
		self.getProperty().SetColor(colorMap.getR(), colorMap.getG(),
			colorMap.getB())


	def setActor(self):
		"""
		Set up the 2D text actor, its mapper and its display position.
		"""

		self.vtk_text_actor = vtk.vtkScaledTextActor()
		self.vtk_text_actor.SetMapper(self.vtk_text_mapper)
		self.vtk_text_actor.SetDisplayPosition(10, 10)

	def addActor(self):
		"""
		Add the 2D text actor to the renderer.
		"""

		self.scene.getRenderer().AddActor2D(self.vtk_text_actor)
