"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common 

class Map(Common):
	"""
	Class that shows a scalar field by color on the domain surface. 
	"""

	def __init__(self, scene, data_collector, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization	
		@type lut: L{BlueToRed <colormap.BlueToRed>} object or 
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		Common.__init__(self, scene, data_collector)
		Common.setMapperInput(self, self.data_collector.getReader().GetOutput(),
			lut)
		Common.setActorInput(self)
		Common.addActor(self)


from plane import Plane 

class MapOnPlane(Plane):
	"""
	Class that shows a scalar field by color on a given plane.
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
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		# "Cut" is used to distinguish cutting from clipping.
		Plane.__init__(self, scene, data_collector,
			data_collector.getReader().GetOutput(), transform, lut, "Cut")

from plane import Plane

class MapOnClip(Plane):
	"""
	Class that shows a scalar field by color on a given clip.
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
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		# "Clip is used to distinguish clipping from cutting.
		Plane.__init__(self, scene, data_collector, 
			data_collector.getReader().GetOutput(), transform, lut, "Clip")

