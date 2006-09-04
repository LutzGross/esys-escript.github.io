"""
Class that shows a scalar field by a surface map. 
"""

import vtk
from common import * 

class Map(Common):
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
		Common.setMapper(self, "self.data_collector.getReader().GetOutput()")
		Common.setActor(self)
		Common.addActor(self)

"""
class MapOnPlane():
shows scalar data by color on a given plane
"""
pass
