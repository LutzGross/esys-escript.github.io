"""
@author: John NGUI
"""

import vtk
from constant import WarpMode

class Warp:
	"""
	Class that defines the deformation of a scalar field.
	"""

	def __init__(self, object, warp_mode):
		"""
		Initialise the warp scalar/vector.

		@type object: vtkPolyData, etc
		@param object: Input for the warp scalar
		@type warp_mode: L{WarpMode <constant.WarpMode>} constant
		@param warp_mode: Mode in which to deform the data
		"""

		self.__object = object

		if(warp_mode == WarpMode.SCALAR): # Deform data with scalar data.
			self.__vtk_warp = vtk.vtkWarpScalar()
		elif(warp_mode == WarpMode.VECTOR): # Deform data with vector data.
			self.__vtk_warp = vtk.vtkWarpVector()

		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the warp scalar/vector.
		"""

		self.__vtk_warp.SetInput(self.__object)

	def setScaleFactor(self, scale_factor):
		"""
		Set the displacement scale factor.

		@type scale_factor: Number
		@param scale_factor: Scale factor of the displacement
		"""

		self.__vtk_warp.SetScaleFactor(scale_factor)

	def _getOutput(self):
		"""
		Return the output of the deformed data. 
		
		@rtype: vtkPointSet
		@return: PointSet data
		"""
		
		return self.__vtk_warp.GetOutput()

