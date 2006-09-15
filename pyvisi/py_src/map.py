"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common 

class Map(Common):
	"""
	Class that shows a scalar field by a surface map. 
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
		@param lut: Color lookup tabl to be used by the mapper
		"""

		Common.__init__(self, scene, data_collector)
		Common.setMapperInput(self, self.data_collector.getReader().GetOutput(),
			lut)
		Common.setActorInput(self)
		Common.addActor(self)


from plane import Plane 

class MapOnPlane(Plane):
	"""
	Class that shows a scalar field on a given plane.
	"""

	def __init__(self, scene, data_collector, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualzation
		"""

		Plane.__init__(self, scene, data_collector,
			data_collector.getReader().GetOutput(), lut)

