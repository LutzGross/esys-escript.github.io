"""
@author: John NGUI
"""

import vtk

class Glyph3D:
	"""
	Class that defines 3D glyph.
	"""

	def __init__(self, object, source, range):
		"""
		Initialise the 3D glyph.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the 3D glyph 
		@type source: vtkPolyData 	
		@param source: Source for the 3D glyph (i.e. Arrow2D, Arrow3D, etc)
		@type range: Two column list
		@param range: Range to map scalar values into 
		"""

		self.__object = object
		self.__source = source
		self.__range = range
		self.__vtk_glyph3D = vtk.vtkGlyph3D()

		self.__setInput()
		self.__setSource()
		self.setScaleModeByVector()	
		#self.setColorModeByVector()
		self.__setVectorModeByVector()
		self.__setClampingOn()
		self.__setScalingOn()
		self.__setOrientOn()
		self.__setRange()
		self.__output = self.__vtk_glyph3D.GetOutput()

	def __setInput(self):
		"""
		Set the input for the 3D glyph.	
		"""

		self.__vtk_glyph3D.SetInput(self.__object)		

	def __setSource(self):
		"""
		Set the source for the 3D glyph.	
		"""

		self.__vtk_glyph3D.SetSource(self.__source)

	def setScaleModeByVector(self):
		"""
		Set the 3D glyph to scale according to the vector.
		"""

		self.__vtk_glyph3D.SetScaleModeToScaleByVector()

	def setScaleModeByScalar(self):
		"""
		Set the 3D glyph to scale according to the scalar.
		"""

		self.__vtk_glyph3D.SetScaleModeToScaleByScalar()

	def _setColorModeByVector(self):
		"""
		Set the 3D glyph to color according to the vector.
		"""

		self.__vtk_glyph3D.SetColorModeToColorByVector()

	def _setColorModeByScalar(self):
		"""
		Set the 3D glyph to color according to the scalar.
		"""

		self.__vtk_glyph3D.SetColorModeToColorByScalar()

	def __setVectorModeByVector(self):
		"""
		Set the 3D glyph vector mode according to the vector.
		"""

		self.__vtk_glyph3D.SetVectorModeToUseVector()

	def setScaleFactor(self, scale_factor):
		"""
		Set the 3D glyph scaling factor.
		
		@type scale_factor: Number
		@param scale_factor: Scale factor
		"""

		self.__vtk_glyph3D.SetScaleFactor(scale_factor)

	def __setClampingOn(self):
		"""
		Enable clamping of "scalar" values to range.
		"""

		self.__vtk_glyph3D.SetClamping(1)

	def __setScalingOn(self):
		"""
		Enable the scaling of the rendered object.	
		"""

		self.__vtk_glyph3D.ScalingOn()

	def __setOrientOn(self):
		"""
		Enable the orienting of the rendered object along vector/normal.	
		"""

		self.__vtk_glyph3D.OrientOn()
		
	def __setRange(self):
		"""
		Set the 3D glyph range to maps scalar values into.
		"""

		self.__vtk_glyph3D.SetRange(self.__range)

	def _getOutput(self):
		"""
		Return the 3D glyph.

		@rtype: vtkPolyData
		@return Polygonal data
		"""

		return self.__output


class TensorGlyph:

	def __init__(self, object, source):
		self.__object = object
		self.__source = source
		self.__vtk_tensor_glyph = vtk.vtkTensorGlyph()

		self.__setInput()
		self.__setSource()

	def __setInput(self):
		"""
		Set the input for the tensor glyph.
		"""

		self.__vtk_tensor_glyph.SetInput(self.__object)

	def __setSource(self):
		"""
		Set the source for the tensor glyph.
		"""
													
		self.__vtk_tensor_glyph.SetSource(self.__source)

	def setScaleFactor(self, scale_factor):
		"""
		Set the tensor glyph scaling factor.
		
		@type scale_factor: Number
		@param scale_factor: Scale factor
		"""

		self.__vtk_tensor_glyph.SetScaleFactor(scale_factor)

	def setMaxScaleFactor(self, max_scale_factor):
		"""
		Set the maximum allowable scale factor.	

		@type max_scale_factor: Number
		@param max_scale_factor: Maximum allowable scale factor.
		"""

		self.__vtk_tensor_glyph.SetMaxScaleFactor(scale_factor)

	def _getOutput(self):
		"""
		Return the tensor glyph.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_tensor_glyph.GetOutput()
