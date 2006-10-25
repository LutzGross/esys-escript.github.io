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
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBLue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		Common.__init__(self, scene, data_collector)
		self.vtk_contour = vtk.vtkContourFilter() 
		self.setContour()

		Common.setMapperInput(self, self.vtk_contour.GetOutput(), lut)
		Common.setActorInput(self)
		Common.addActor(self)

	def setContour(self):
		"""
		Set up the contour.
		"""

		self.vtk_contour.SetInput(self.data_collector.getReader().GetOutput())

	def generateValues(self, number_contours, min_range, max_range):
		"""
		Generate the specified number of contours within the specified range.

		@type number_contours: Number
		@param number_contours: Number of contours to generate	
		@type min_range: Number
		@param min_range: Minimum contour range 
		@type max_range: Number 
		@param max_range: Maximum contour range 
		"""

		self.vtk_contour.GenerateValues(number_contours, min_range, max_range)

from contour import Contour
from plane import Plane

class ContourOnPlane(Contour, Plane):
	"""
	Class that shows a scalar field by contour surfaces on a given plane.
	"""

	def __init__(self, scene, data_collector, transform, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object		
		@param scene: Scene in which components rae to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@param transform: L{Transform <geo.Transform>} object
		@type transform: Orientation of the plane
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBLue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		# Declared because they are needed by the setContour method.
		self.data_collector = data_collector
		self.vtk_contour = vtk.vtkContourFilter() 
	
		Contour.setContour(self)
		# "Cut" is used to distinguish cutting from clipping.
		Plane.__init__(self, scene, data_collector,
			self.vtk_contour.GetOutput(), transform, lut, "Cut")


from contour import Contour
from plane import Plane

class ContourOnClip(Contour, Plane):
	"""	
	Class that shows a scalar field by contour surfaces on a given clip.
	"""

	def __init__(self, scene, data_collector, transform, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object		
		@param scene: Scene in which components rae to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@param transform: L{Transform <geo.Transform>} object
		@type transform: Orientation of the plane
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBLue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		# Declared because they are needed by the setContour method.
		self.data_collector = data_collector
		self.vtk_contour = vtk.vtkContourFilter() 

		 
		Contour.setContour(self)
		# "Clip" is used to distinguish clipping from cutting.
		Plane.__init__(self, scene, data_collector,
			self.vtk_contour.GetOutput(), transform, lut, "Clip")
