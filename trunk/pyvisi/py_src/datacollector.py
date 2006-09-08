"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import *
from colormap import ColorMap

class DataCollector(Common):
	"""
	Class that deals with data for the visualization.
	"""

	def __init__(self, scene, outline = True):
		"""
		Initialize all the instance variables. 

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type outline: Boolean (I{True or False})
		@param outline: Determines the outline for the rendered object
		"""

		self.scene = scene	
		self.outline = True
		self.file_name = None
		self.vtk_outline = None
		self.vtk_xml_reader = None 
		self.vtk_xml_reader_output = None

	def setFileName(self, file_name):
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


