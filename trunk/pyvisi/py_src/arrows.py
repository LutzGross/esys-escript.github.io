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
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
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
		# Default vector mode is vector.
		self.setVectorMode("Vector")
		# Default scale mode is vector.
		self.setScaleMode("Vector")
		# Default color mode is scalar.
		self.setColorMode("Scalar")
		self.setScaleFactor(0.2) 

	def setVectorMode(self, vector_mode):
		"""
		Set the arrows vector mode.
		@type vector_mode: String
		@param vector_mode: Arrows vector mode
		"""

		eval("self.vtk_glyph.SetVectorModeToUse%s" % vector_mode)

	def setScaleMode(self, scale_mode):
		"""
		Set the arrows scale mode.
		@type scale_mode: String
		@param scale_mode: Arrows scale mode
		"""

		eval("self.vtk_glyph.SetScaleModeToScaleBy%s" % scale_mode)

	def setScaleFactor(self, scale_factor):
		"""
		Set the arrows scale factor.
		@type scale_factor: Number
		@param scale_factor: Size of the arrows 
		"""

		self.vtk_glyph.SetScaleFactor(scale_factor)

	def setColorMode(self, color_mode):
		"""
		Set the arrows color mode.	
		@type color_mode: String
		@param color_mode: Arrows color mode
		"""
	
		eval("self.vtk_glyph.SetColorModeToColorBy%s()" % color_mode)


from arrows import Arrows
from plane import Plane 
	
class ArrowsOnPlane(Arrows, Plane):
	"""
	Class that shows a vector field by arrows on a given plane.
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

		# Declared because they are needed by the setArrows method.
		self.data_collector = data_collector
		self.vtk_glyph = vtk.vtkGlyph3D()
		
		Arrows.setArrows(self)
		# "Cut" is used to distinguish cutting from clipping.
		Plane.__init__(self, scene, data_collector,
			self.vtk_glyph.GetOutput(), transform, lut, "Cut")

from arrows import Arrows
from plane import Plane

class ArrowsOnClip(Arrows, Plane):
	"""
	Class that shows a vector field by arrows on a clip.
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

		# Declared because they are needed by the setArrows method.
		self.data_collector = data_collector
		self.vtk_glyph = vtk.vtkGlyph3D()
		
		Arrows.setArrows(self)
		# "Clip" is used to distinguish clipping from cutting.
		Plane.__init__(self, scene, data_collector,
			self.vtk_glyph.GetOutput(), transform, lut, "Clip")

