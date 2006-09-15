"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common 
from constants import *
from style import Style

class DataCollector(Common):
	"""
	Class that deals with data for the visualization.
	"""

	def __init__(self, scene, outline = True, cube_axes = True):
		"""
		Initialize all the instance variables. 

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type outline: Boolean (I{True or False})
		@param outline: Determines the outline for the rendered object
		"""

		Common.__init__(self, scene)
		self.outline = outline 
		self.cube_axes = cube_axes
		self.vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()

	def setFileName(self, file_name):
		"""	
		Set up the file reader and set the file name.

		@type file_name: String
		@param file_name: Name of the file to be read.
		"""

		self.vtk_xml_reader.SetFileName(file_name)

		if(self.outline == True):
			self.setOutline()

		if(self.cube_axes == True):
			self.setCubeAxes()

	def setOutline(self):
		"""	
		Set the outline for the rendered object.
		"""

		self.vtk_outline = vtk.vtkOutlineFilter()
		self.vtk_outline.SetInput(self.vtk_xml_reader.GetOutput())	
	
		Common.setMapperInput(self, self.vtk_outline.GetOutput())
		Common.setActorInput(self)
		Common.addActor(self)
		# Default outline is black color.	
		Common.setActorColor(self, BLACK) 

	def setCubeAxes(self):
		vtk_cube_axes = vtk.vtkCubeAxesActor2D()
		vtk_cube_axes.SetInput(self.getReader().GetOutput())
		vtk_cube_axes.SetCamera(self.scene.getRenderer().GetActiveCamera())
		vtk_cube_axes.SetLabelFormat("%6.4g")
		vtk_cube_axes.SetFlyModeToOuterEdges()
		vtk_cube_axes.SetFontFactor(0.9)

		style = Style() 		
		style.setColor(BLACK)
		style.setShadow()

		vtk_cube_axes.SetAxisTitleTextProperty(style.getStyle())
		vtk_cube_axes.SetAxisLabelTextProperty(style.getStyle())

		self.scene.getRenderer().AddActor(vtk_cube_axes)

	def getReader(self):
		"""
		Return the file reader.

		@rtype: vtkXMLUnstructuredGridReader
		@return: VTK XML unstructured grid reader
		"""

		return self.vtk_xml_reader

