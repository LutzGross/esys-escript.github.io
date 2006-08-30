"""

classes dealing with data for the visualization

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane, L. Gross"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision:$"
__date__="$Date:$"

import vtk
from common import *

class DataCollector(Common):

	def __init__(self, open_scene, outline = True):
		self.open_scene = open_scene	
		self.outline = True
		self.file_name = None
		self.vtk_outline = None
		self.vtk_xml_reader = None 
		self.vtk_xml_reader_output = None

	# set up the file reader and set the file name
	def setSource(self, file_name):
		self.file_name = file_name
		self.vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
		self.vtk_xml_reader.SetFileName(self.file_name)

		if(self.outline == True):
			self.setOutline()
			Common.setMapper(self, "self.vtk_outline.GetOutput()")
			Common.setActor(self)
			Common.addActor(self)	
			Common.setColor(self, 0, 0, 0)
			
			
	
	# return the file reader output
	def getReader(self):
		return self.vtk_xml_reader

	# set the outline
	def setOutline(self):
		self.vtk_outline = vtk.vtkOutlineFilter()
		self.vtk_outline.SetInput(self.vtk_xml_reader.GetOutput())	
		
		#vtk_outline_mapper = vtk.vtkPolyDataMapper()
		#vtk_outline_mapper.SetInput(vtk_outline.GetOutput())
		
		#vtk_outline_actor = vtk.vtkActor()
		#vtk_outline_actor.SetMapper(vtk_outline_mapper)
		#vtk_outline_actor.GetProperty().SetColor(0, 0, 0)	

		#self.open_scene.getRenderer().AddActor(vtk_outline_actor)	
	

		
		
						


