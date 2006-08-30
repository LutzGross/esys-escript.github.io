"""
class that shows scalar data by contour surfaces

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

class Contour:
	
	def __init__(self, open_scene, data_collector):
		self.open_scene = open_scene
		self.data_collector = data_collector
		self.vtk_contour = None
		self.vtk_contour_mapper = None
		self.vtk_contour_actor = None

		self.setContour()
		self.setMapper()
		self.setActor()

	# set up the contour and specify the number and range
	def setContour(self):
		self.vtk_contour = vtk.vtkContourFilter()
		self.vtk_contour.SetInput(self.data_collector.getReader().GetOutput())
		self.vtk_contour.GenerateValues(5, 0.0, 1.2)
		
	# set up the mapper and data	
	def setMapper(self):
		self.vtk_contour_mapper = vtk.vtkPolyDataMapper()
		self.vtk_contour_mapper.SetInput(
			self.vtk_contour.GetOutput())
	
	# set up the actor and add the actor to the scene
	def setActor(self):
		self.vtk_contour_actor = vtk.vtkActor()
		self.vtk_contour_actor.SetMapper(self.vtk_contour_mapper)
		self.vtk_contour_actor.GetProperty().SetOpacity(0.6)

		self.open_scene.getRenderer().AddActor(self.vtk_contour_actor)


#class ContourOnPlane(Component):
"""
shows scalar data by contour surfaces on a given plane
"""
pass
