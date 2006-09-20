"""
@author: John Ngui
@author: Lutz Gross	
"""

import vtk
from common import *

class Contour(Common):
	"""	
	Class that shows a scalar field by contour surfaces.
	"""

	def __init__(self, scene, data_collector, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		"""

		Common.__init__(self, scene, data_collector)
		self.vtk_contour = vtk.vtkContourFilter() 
		self.setContour()

		Common.setMapperInput(self, self.vtk_contour.GetOutput(), lut)
		Common.setActorInput(self)
		Common.addActor(self)

	def setContour(self):
		"""
		Set up the contour and its input.
		"""

		self.vtk_contour.SetInput(self.data_collector.getReader().GetOutput())


	def generateValues(self, number_contours, min_range, max_range):
		"""
		Generate the specified number of contours within the specified range.

		@type number_contours: Number
		@param number_contours: Number of contours to generate	
		@type min_range: Number
		@param min_range: Minimum contour value
		@type max_range: Number 
		@param max_range: Maximum contour value
		"""

		self.vtk_contour.GenerateValues(number_contours, min_range, max_range)


from contour import Contour
from plane import Plane

class ContourOnPlane(Contour, Plane):
	"""
	Class that shows a scalar field by contour surfaces on a given plane.
	"""

	def __init__(self, scene, data_collector, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object		
		@param scene: Scene in which components rae to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		"""

		self.data_collector = data_collector
		self.vtk_contour = vtk.vtkContourFilter() 
		Contour.setContour(self)

		Plane.__init__(self, scene, data_collector,
			self.vtk_contour.GetOutput(), lut)


from contour import Contour 

class IsoSurface(Contour):

	def __init__(self, scene, data_collector, lut = None):
		Contour.__init__(self, scene, data_collector, lut)

	def setValue(self, contour_number, value):
		"""
		Set the contour number and its value.

		@type contour_number: Number
		@param contour_number: Contour number
		@type value: Number
		@param value: Contour value
		"""

		self.vtk_contour.SetValue(contour_number, value)

from contour import IsoSurface, ContourOnPlane

class IsoSurfaceOnPlane(IsoSurface, ContourOnPlane):

	def __init__(self, scene, data_collector, lut = None):
		ContourOnPlane.__init__(self, scene, data_collector, lut)
		
