"""
@author: John NGUI
"""

import vtk

class Clipper:
	"""
	Class that defines a clipper.
	"""
	
	# plane by default is assigned to None because a plane is not required 
	# when a scalar value is used to perform the clipping.
	def __init__(self, object, plane = None):
		"""
		Initialise the clipper.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the clipper
		@type plane: vtkPlane
		@param plane: Plane to clip the object
		"""

		self.__object = object 
		if(plane != None): # True only if a plane is used to perform clipping.
			self.__plane = plane
			
		self.__vtk_clipper = vtk.vtkClipDataSet()

		self.__setInput()
		#self.__setClipFunction()
		# Generates the clipped away section. At this stage, the clipped
		# away section is not used.
		#self.__vtk_clipper.GenerateClippedOutputOn()
		self.setInsideOutOn()
		self.__vtk_clipper.Update()

	def __setInput(self):
		"""
		Set the input for the clipper.
		"""

		self.__vtk_clipper.SetInput(self.__object)

	def _setClipFunction(self):
		"""
		Set the clip function.
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

