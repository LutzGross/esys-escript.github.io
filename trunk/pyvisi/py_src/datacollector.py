"""
Class that deals with data for the visualization.
"""

import vtk
from common import *

class DataCollector(Common):
	"""
	@author: John Ngui
	@author: Lutz Gross
	"""

	def __init__(self, open_scene, outline = True):
		"""
		Initialize all the instance variables. 

		@type open_scene: L{OpenScene <openscene.OpenScene>} object
		@param open_scene: Scene in which components are to be added to
		@type outline: Boolean (I{True or False})
		@param outline: Determines the outline for the rendered object
		"""

		self.open_scene = open_scene	
		self.outline = True
		self.file_name = None
		self.vtk_outline = None
		self.vtk_xml_reader = None 
		self.vtk_xml_reader_output = None

	def setSource(self, file_name):
		"""	
		Set up the file reader and set the file name.

		@type file_name: String
		@param file_name: Name of the file to be read.
		"""

		self.file_name = file_name
		self.vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
		self.vtk_xml_reader.SetFileName(self.file_name)

		if(self.outline == True):
			self.setOutline()
			Common.setMapper(self, "self.vtk_outline.GetOutput()")
			Common.setActor(self)
			Common.addActor(self)	
			Common.setColor(self, 0, 0, 0) # Default outline is black
	
	def getReader(self):
		"""
		Return the file reader.

		@rtype: vtkXMLUnstructuredGridReader
		@return: VTK XML unstructured grid reader
		"""

		return self.vtk_xml_reader

	def setOutline(self):
		"""	
		Set the outline for the rendered object.
		"""

		self.vtk_outline = vtk.vtkOutlineFilter()
		self.vtk_outline.SetInput(self.vtk_xml_reader.GetOutput())	


