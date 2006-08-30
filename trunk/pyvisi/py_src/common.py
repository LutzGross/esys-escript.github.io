"""
class that shows a vector field by arrows

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

class Common:

	def __init__(self, open_scene, data_collector):
		self.open_scene = open_scene
		self.data_collector = data_collector
		self.vtk_mapper = None
		self.vtk_actor = None

	def setMapper(self, component):
		self.vtk_mapper = vtk.vtkDataSetMapper()
		eval("self.vtk_mapper.SetInput(%s)" % component)

	def setActor(self):
		self.vtk_actor = vtk.vtkActor()
		self.vtk_actor.SetMapper(self.vtk_mapper)

	def addActor(self):
		self.open_scene.getRenderer().AddActor(self.vtk_actor) 

	def setOpacity(self, opacity):
		self.getProperty().SetOpacity(opacity)

	def setColor(self, red, green, blue):
		self.getProperty().SetColor(red, green, blue)

	def setRepresentation(self, representation):
		eval("self.getProperty().SetRepresentationTo%s()" % representation)
	
	def getProperty(self):
		return self.vtk_actor.GetProperty()

		
