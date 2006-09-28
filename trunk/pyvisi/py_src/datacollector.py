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
	Class that deals with data for visualization.
	"""

	def __init__(self, scene, outline = True, cube_axes = True):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type outline: Boolean (I{True or False})
		@param outline: Places or removes an outline for the rendered object
		@type cube_axes: Boolean ({True or False})
		@param cube_axes: Places or removes a cube axes for the rendered object
		"""

		Common.__init__(self, scene)
		self.outline = outline 
		self.cube_axes = cube_axes
		self.vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()

	def setFileName(self, file_name):
		"""	
		Set up the file reader and file name from which data is to be read.
		@type file_name: String
		@param file_name: Name of the file to be read.
		"""

		self.vtk_xml_reader.SetFileName(file_name)

		# Set up the outline and cube axes only if true.
		if(self.outline == True):
			self.setOutline()
		if(self.cube_axes == True):
			self.setCubeAxes()

	def setOutline(self):
		"""	
		Set up the outline for the rendered object.
		"""

		self.vtk_outline = vtk.vtkOutlineFilter()
		self.vtk_outline.SetInput(self.vtk_xml_reader.GetOutput())	
	
		Common.setMapperInput(self, self.vtk_outline.GetOutput())
		Common.setActorInput(self)
		Common.addActor(self)
		# Default outline is black.	
		Common.setActorColor(self, BLACK) 

	def setCubeAxes(self):
		"""	
		Set up the cube axes for the rendered object.
		"""

		vtk_cube_axes = vtk.vtkCubeAxesActor2D()
		vtk_cube_axes.SetInput(self.getReader().GetOutput())
		vtk_cube_axes.SetCamera(self.scene.getRenderer().GetActiveCamera())
		# Formats the labes on the axes.
		vtk_cube_axes.SetLabelFormat("%6.4g")
		# Use the outer edges of the bounding box to draw the axes.
		vtk_cube_axes.SetFlyModeToOuterEdges()
		vtk_cube_axes.SetFontFactor(0.9)

		# Style for the axes title and label.
		style = Style() 		
		style.setColor(BLACK)
		style.shadowOn()

		vtk_cube_axes.SetAxisTitleTextProperty(style.getStyle())
		vtk_cube_axes.SetAxisLabelTextProperty(style.getStyle())

		self.scene.getRenderer().AddActor(vtk_cube_axes)

	def getReader(self):
		"""
		Return the file reader.

		@rtype: vtkXMLUnstructuredGridReader
		@return: VTK XML unstructured grid reader that is used to read files.
		"""

		return self.vtk_xml_reader

