"""
class that represents scalar data as plane deformated along the plane normal 
and proportional to the scalar value on the plane.
"""

import vtk
from plane import Plane
from common import Common 

class Carpet(Common, Plane):
	"""
	represents scalar data as plane deformated along the plane normal 
	and proportional to the scalar value on the plane.
	"""

	def __init__(self, scene, data_collector, lut = None):
		Common.__init__(self, scene, data_collector)
		self.vtk_plane = vtk.vtkPlane()
		self.vtk_cutter = vtk.vtkCutter()
		self.vtk_transform = vtk.vtkTransform()
		self.vtk_transform_filter = vtk.vtkTransformPolyDataFilter()
		self.vtk_warp = vtk.vtkWarpScalar()
		
		Plane.setPlane(self)
		Plane.setCutter(self, data_collector.getReader().GetOutput())	
		self.warpScalar()
		Plane.setTransformFilter(self, self.vtk_warp.GetOutput())	
		
		Common.setMapperInput(self, self.vtk_transform_filter.GetOutput(), lut)
		Common.setActorInput(self)
		Common.addActor(self)

	def warpScalar(self):
		self.vtk_warp.SetInput(self.vtk_cutter.GetOutput())
		self.vtk_warp.SetScaleFactor(0.5)
			
							


		
	
