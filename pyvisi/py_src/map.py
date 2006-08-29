"""
class that shows scalar data by color on the domain surface

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

class Map:

	def __init__(self, scene, data_collector):
		self.scene = scene
		self.data_collector = data_collector
		self.vtk_xml_mapper = None
		self.vtk_xml_actor = None

		self.setMapper()
		self.setActor()

	def setMapper(self):
		self.vtk_xml_mapper = vtk.vtkDataSetMapper()
		self.vtk_xml_mapper.SetInput(
			self.data_collector.getReader().GetOutput())
	
	def setActor(self):
		self.vtk_xml_actor = vtk.vtkActor()
		self.vtk_xml_actor.SetMapper(self.vtk_xml_mapper)

		self.scene.getRenderer().AddActor(self.vtk_xml_actor)

	

"""
class MapOnPlane():
shows scalar data by color on a given plane
"""
pass
