"""
Class that shows a scalar field by contour surfaces.
"""

import vtk
from common import *

class Contour(Common):
	"""
	@author: John Ngui
	@author: Lutz Gross	
	"""
	
	def __init__(self, open_scene, data_collector):
		"""
		@type open_scene: L{OpenScene <openscene.OpenScene>} object
		@param open_scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		"""

		Common.__init__(self, open_scene, data_collector)
		self.setContour()

		Common.setMapper(self, "self.vtk_contour.GetOutput()")
		Common.setActor(self)
		Common.addActor(self)

	def setContour(self):
		"""
		Set up the contour and its input.
		"""

		self.vtk_contour = vtk.vtkContourFilter()
		self.vtk_contour.SetInput(self.data_collector.getReader().GetOutput())

	def setValue(self, contour_number, value):
		"""
		Set the contour number and its value.

		@type contour_number: Number
		@param contour_number: Contour number
		@type value: Number
		@param value: Contour value
		"""

		self.vtk_contour.SetValue(contour_number, value)

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


#class ContourOnPlane(Component):
"""
shows scalar data by contour surfaces on a given plane
"""
pass
