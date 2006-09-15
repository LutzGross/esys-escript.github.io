"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common 

class Arrows(Common):
	"""
	Class that shows a vector field by arrows.
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
		self.vtk_glyph = vtk.vtkGlyph3D()
		self.setArrows()

		Common.setMapperInput(self, self.vtk_glyph.GetOutput(), lut)
		Common.setActorInput(self)
		Common.addActor(self)		
	
	def setArrows(self):
		"""
		Set up the glyph and use arrows as the source.
		"""

		vtk_arrows = vtk.vtkArrowSource()
		
		self.vtk_glyph.SetInput(self.data_collector.getReader().GetOutput())
		self.vtk_glyph.SetSource(vtk_arrows.GetOutput())
		self.setVectorMode("Vector")
		self.setScaleMode("Vector")
		self.setColorMode("Scalar")
		self.setScaleFactor(0.3) 

	def setVectorMode(self, vector_mode):
		eval("self.vtk_glyph.SetVectorModeToUse%s" % vector_mode)

	def setScaleMode(self, scale_mode):
		eval("self.vtk_glyph.SetScaleModeToScaleBy%s" % scale_mode)

	def setScaleFactor(self, scale_factor):
		"""
		Set the scale factor for the arrows.

		@type scale_factor: Number
		@param scale_factor: Scale factor
		"""

		self.vtk_glyph.SetScaleFactor(scale_factor)

	def setColorMode(self, color_mode):
		"""
		Set the color mode for the arrows.	

		@type color_mode: String
		@param color_mode: Color mode for the arrows (I{Scalar or Vector})
		"""
	
		eval("self.vtk_glyph.SetColorModeToColorBy%s()" % color_mode)


from arrows import Arrows
from plane import Plane 
	
class ArrowsOnPlane(Arrows, Plane):
	"""
	Class that shows a vector field by arrows on a given plane.
	"""

	def __init__(self, scene, data_collector, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>} 
			object

		@param data_collector: Source of data for visualization
		"""

		self.data_collector = data_collector
		self.vtk_glyph = vtk.vtkGlyph3D()
		Arrows.setArrows(self)

		Plane.__init__(self, scene, data_collector,
			self.vtk_glyph.GetOutput(), lut)

