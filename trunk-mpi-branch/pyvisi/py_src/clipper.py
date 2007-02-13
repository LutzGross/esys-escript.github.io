"""
@author: John NGUI
"""

import vtk

class Clipper:
	"""
	Class that defines a clipper.
	"""
	
	def __init__(self, object, plane):
		"""
		Initialise the clipper.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the clipper
		@type plane: vtkPlane
		@param plane: Plane to clip the rendered object
		"""

		self.__object = object 
		# True only if a plane is used to perform clipping. Will be false 
		# for scalar clipping.
		if(plane != None): 
			self.__plane = plane
			
		self.__vtk_clipper = vtk.vtkClipDataSet()
		self.__setupClipper()

	def __setupClipper(self):
		"""
		Setup the clipper.
		"""

		self.__setInput()
		self.setInsideOutOn()
		self.__vtk_clipper.Update()

	def __setInput(self):
		"""
		Set the input for the clipper.
		"""

		self.__vtk_clipper.SetInput(self.__object)

	def _setClipFunction(self):
		"""
		Set the clip function (using a plane).
		"""

		self.__vtk_clipper.SetClipFunction(self.__plane)

	def setInsideOutOn(self):
		"""
		Clips the one side of the rendered object.		
		"""

		self.__vtk_clipper.InsideOutOn()

	def setInsideOutOff(self):
		"""
		Clips the other side of the rendered object.		
		"""

		self.__vtk_clipper.InsideOutOff()

	def setClipValue(self, value):
		"""
		Set the scalar clip value.

		@type value: Number
		@param value: Scalar clip value
		"""

		self.__vtk_clipper.SetValue(value)

	def _getOutput(self):
		"""
		Return the output of the clipper.

		@rtype: vtkUnstructuredGrid
		@return: Unstructured grid
		"""

		return self.__vtk_clipper.GetOutput()

