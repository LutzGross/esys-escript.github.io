"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common

class Tensor(Common):
	"""
	Class that show a tensor field by ellipsoids.
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
		self.vtk_tensor_glyph = vtk.vtkTensorGlyph()
		self.vtk_sphere = vtk.vtkSphereSource()
		self.vtk_poly_data_normals = vtk.vtkPolyDataNormals()
		self.setTensor()

		Common.setMapperInput(self, self.vtk_poly_data_normals.GetOutput(), lut)
		Common.setActorInput(self)
		Common.addActor(self)

	def setTensor(self):
		"""
		Set up the tensor glyph and use sphere as the source.		
		"""

		# Default theta and phi resolution is 10.
		self.setThetaResolution(10)
		self.setPhiResolution(10)
	
		self.vtk_tensor_glyph.SetInput(
			self.data_collector.getReader().GetOutput())
		# Use sphere are the source for the glyph.
		self.vtk_tensor_glyph.SetSource(self.vtk_sphere.GetOutput())
		# Default scale and max scale factor is 0.2 and 5 respectively.
		self.setScaleFactor(0.2)
		self.setMaxScaleFactor(5)

		# Needed to generate the right colors. If not used, some glyphs
		# will be black in color.
		self.vtk_poly_data_normals.SetInput(self.vtk_tensor_glyph.GetOutput())
		
	def setThetaResolution(self, resolution):
		"""
		Set the number of points in the longitude direction.
		@type resolution: Number
		@param resolution: Number of points in the longitude direction
		"""

		self.vtk_sphere.SetThetaResolution(resolution)

	def setPhiResolution(self, resolution):
		"""
		Set the number of points in the latitude direction.
		@type resolution: Number
		@param resolution: Number of points in the latitude direction
		"""

		self.vtk_sphere.SetPhiResolution(resolution)

	def setScaleFactor(self, scale_factor):
		"""
		Set the tensor glyph scale factor.
		@type scale_factor: Number
		@param scale_factor: Size of the sphere
		"""

		self.vtk_tensor_glyph.SetScaleFactor(scale_factor)

	def setMaxScaleFactor(self, max_scale_factor):
		"""
		Set the maximum allowable scale factor.
		@type max_scale_factor: Number
		@param max_scale_factor: Maximum size of the sphere
		"""

		self.vtk_tensor_glyph.SetMaxScaleFactor(max_scale_factor)

from tensor import Tensor
from plane import Plane

class TensorOnPlane(Tensor, Plane):
	"""
	Class that show a tensor field by ellipsoids on given a plane.
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

		# Declared because they are needed by the setTensor method.
		self.data_collector = data_collector
		self.vtk_tensor_glyph = vtk.vtkTensorGlyph()
		self.vtk_sphere = vtk.vtkSphereSource()
		self.vtk_poly_data_normals = vtk.vtkPolyDataNormals()

		Tensor.setTensor(self)
		# "Cut" is used to distinguish cutting from clipping.
		Plane.__init__(self, scene, data_collector, 
			self.vtk_poly_data_normals.GetOutput(), transform, lut, "Cut")

from tensor import Tensor
from plane import Plane

class TensorOnClip(Tensor, Plane):
	"""
	Class that show a tensor field by ellipsoids on a given clip.
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
		
		# Declared because they are needed by the setTensor method.
		self.data_collector = data_collector
		self.vtk_tensor_glyph = vtk.vtkTensorGlyph()
		self.vtk_sphere = vtk.vtkSphereSource()
		self.vtk_poly_data_normals = vtk.vtkPolyDataNormals()
		
		Tensor.setTensor(self)
		# "Clip" is used to distinguish clipping from cutting.
		Plane.__init__(self, scene, data_collector, 
			self.vtk_poly_data_normals.GetOutput(), transform, lut, "Clip")
	

