import vtk
from plane import Plane
from common import Common 

class Carpet(Common, Plane):
	"""
	Class that represents a scalar field as a plane deformated along the plane 
	normal and proportional to the scalar value on the plane.
	"""

	def __init__(self, scene, data_collector, transform, lut = None, 
		deform = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@type transform: L{Transform <geo.Transform>} object
		@param transform: Orientation of the plane
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		@type deform: String 
		@param deform: Mode the data is deformed. Either by I{Scalar} 
			or I{Vector}
		"""

		Common.__init__(self, scene, data_collector)
		# Declared because needed by the setPlane method.
		self.vtk_plane = vtk.vtkPlane()
		self.vtk_cutter = vtk.vtkCutter()
		self.transform = transform.getTransform() 
		
		if(deform == "Scalar"):
			self.vtk_warp = vtk.vtkWarpScalar()
		else:
			self.vtk_warp = vtk.vtkWarpVector()
			
		
		Plane.setPlane(self)
		Plane.setCutter(self, data_collector.getReader().GetOutput())	
		self.warpScalar()
		
		#Common.setMapperInput(self, self.vtk_warp.GetOutput(), lut)
		Common.setMapperInput(self, self.vtk_cutter.GetOutput(), lut)
		Common.setActorInput(self)
		Common.addActor(self)

	def warpScalar(self):
		"""
		Set up the warp scalar and deform the plane with scalar data.
		"""

		self.vtk_warp.SetInput(self.vtk_cutter.GetOutput())

	def setScaleFactor(self, scale_factor):
		"""
		Set the displacement scale factor.
		@type scale_factor: Number
		@param scale_factor: Size of the displacement
		"""

		self.vtk_warp.SetScaleFactor(scale_factor)
