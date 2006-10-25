"""
@author: John Ngui
@author: Lutz Gross
"""

from contour import Contour

class IsoSurface(Contour):
	"""
	Class that shows a scalar field for a given value by an isosurface.
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
		
		Contour.__init__(self, scene, data_collector, lut)

	def setValue(self, contour_number, value):
		"""
		Set the contour number and value.
		@type contour_number: Number
		@param contour_number: Contour number
		@type value: Number
		@param value: Contour value
		"""
		
		self.vtk_contour.SetValue(contour_number, value)

from isosurface import IsoSurface
from contour import  ContourOnPlane

class IsoSurfaceOnPlane(IsoSurface, ContourOnPlane):
	"""
	Class that shows a scalar field for a given value by an isosurface
	on a given plane.
	"""

	def __init__(self, scene, data_collector, transform, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@type transform: L{Transform <geo.Transform>} object
		@param transform: Orientation of the plane
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBLue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		ContourOnPlane.__init__(self, scene, data_collector, transform, lut)

from isosurface import IsoSurface
from contour import  ContourOnClip

class IsoSurfaceOnClip(IsoSurface, ContourOnClip):
	"""
	Class that shows a scalar field for a given value by an isosurface
	on a given clip.
	"""

	def __init__(self, scene, data_collector, transform, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@type transform: L{Transform <geo.Transform>} object
		@param transform: Orientation of the plane
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBLue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		ContourOnClip.__init__(self, scene, data_collector, transform, lut)

